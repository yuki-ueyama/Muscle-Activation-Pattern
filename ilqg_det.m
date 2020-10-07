function [x, u, cost] ...
    = ilqg_det(fnDyn, fnCost, dt, n, x0, u0, varargin)
global flgFix;

%---------------------- User-adjustable parameters ------------------------
% Default
lambdaInit = 0.1;        % Initial value of Levenberg-Marquardt lambda
lambdaFactor = 10;       % Factor for multiplying or dividing lambda
% For stabilization control (Case 5) and sensitivity analysis
% lambdaInit = 1;
% lambdaFactor = 5;     
%  For sensitivity analysis
lambdaInit = 1000;
lambdaFactor = 5;

lambdaMax = 1e+10;       % Exit if lambda exceeds threshold
epsConverge = 1e-10;     % Exit if relative improvement below threshold
maxIter = 5000;          % Exit if number of iterations reaches threshold
flgPrint = 1;            % Show cost- 0:never, 1:every iter, 2:final
maxValue = Inf;          % Upper bound on states and costs (Inf: none)
minIter = 100;           % Exit if number of iterations overs threshold
%---------------------------------------------- Get optional arguments
if nargin>=7
    uMin = varargin{1};
else
    uMin = -Inf;
end
if nargin>=8
    uMax = varargin{2};
else
    uMax = Inf;
end

%---------------------------------------------- Initialization
szX = size(x0, 1);          % Size of state vector
szU = size(u0, 1);          % Size of control vector
L = zeros(szU, szX, n-1);   % Init feedback gains

if size(u0,2)==1            % Init control sequence
    u = repmat(u0, [1,n-1]);
else
    u = u0;
end

if isscalar(uMin)           % Transform scalar arguments into vectors
    uMin = repmat(uMin, [szU,1]);
end
if isscalar(uMax) 
    uMax = repmat(uMax, [szU,1]);
end

flgFix = 0;                 % Clear large-fix global flag
          
% init noise paremeters
A = zeros(szX, szX, n-1);
B =zeros(szX, szU, n-1);

% Init state sequence and cost
[x, cost] = simulate(fnDyn, fnCost, dt, x0, u, maxValue);

%---------------------------------------------- optimization loop
lambda = lambdaInit;
flgChange = 1;

s0 = zeros(1, n);
s = zeros(szX, n);
S = zeros(szX, szX, n);

for iter = 1:maxIter
    %------ STEP 1: approximate dynamics and cost along new trajectory
    if flgChange
        [s0(n),s(:,n),S(:,:,n)] = fnCost(x(:,n), NaN, NaN);  % final cost
        
        for k = n-1:-1:1 % start:increment:end
            % quadratize cost, adjust for dt
            [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x(:,k), u(:,k), k);
            q0(k) = dt * l0;
            q(:,k) = dt * l_x;
            Q(:,:,k) = dt * l_xx;
            r(:,k) = dt * l_u;
            R(:,:,k) = dt * l_uu;
            P(:,:,k) = dt * l_ux;            
            
            % linearize dynamics, adjust for dt
            [~, f_x, f_u] = ...
                fnDyn(x(:, k), u(:, k));
            A(:,:,k) = eye(szX) + dt * f_x;
            B(:,:,k) = dt * f_u;
            
        end        
        flgChange = 0;
    end
    
    
    %------ STEP 2: compute optimal control law and cost-to-go
    for k = n-1:-1:1
        % compute shortcuts g,G,H
        g = r(:,k) + B(:,:,k)'*s(:,k+1);
        G = P(:,:,k) + B(:,:,k)'*S(:,:,k+1)*A(:,:,k);
        H = R(:,:,k) + B(:,:,k)'*S(:,:,k+1)*B(:,:,k);
        
        % find control law
        [l(:,k), L(:,:,k)] = uOptimal(g,G,H,u(:,k),uMin,uMax,lambda);
        
        % update cost-to-go approximation
        S(:,:,k) = Q(:,:,k) + A(:,:,k)'*S(:,:,k+1)*A(:,:,k) + ...
            L(:,:,k)'*H*L(:,:,k) + L(:,:,k)'*G + G'*L(:,:,k);
        s(:,k) = q(:,k) + A(:,:,k)'*s(:,k+1) + ...
            L(:,:,k)'*H*l(:,k) + L(:,:,k)'*g + G'*l(:,k);
        s0(k) = q0(k) + s0(k+1) + l(:,k)'*H*l(:,k)/2 + l(:,k)'*g;
        
        % HACK USED TO PREVENT OCCASIONAL DIVERGENCE
        if isfinite(maxValue)
            S(:,:,k) = fixBig(S(:,:,k), maxValue);
            s(:,k) = fixBig(s(:,k), maxValue);
            s0(k) = fixBig(s0(k), maxValue);
        end
    end    
    
    %------ STEP 3: new control sequence, trajectory, cost
    % simulate linearized system to compute new control
    dx = zeros(szX, 1);
    unew = zeros(szU, n-1);
    for k=1:n-1
        du = l(:,k) + L(:,:,k)*dx;
        du = min(max(du+u(:,k),uMin), uMax) - u(:,k);
        dx = A(:,:,k)*dx + B(:,:,k)*du;
        unew(:,k) = u(:,k) + du;
        
    end
    
    % simulate system to compute new trajectory and cost
    [xnew, costnew] = ...
        simulate(fnDyn, fnCost, dt, x0, unew, maxValue);
    
    %------ STEP 4: Levenberg-Marquardt method
    if costnew<cost
        % decrease lambda (get closer to Newton method)
        lambda = lambda / lambdaFactor;
        
        % accept changes, flag changes
        u = unew;
        x = xnew;
        
        flgChange = 1;
        
        if flgPrint==1
            fprintf('Iteration = %d;  Cost = %.6f;  logLambda = %.1f\n', ...
                iter, costnew, log10(lambda) );
        end
        
        if iter>minIter && (abs(costnew - cost)/cost < epsConverge)
            cost = costnew;
            break;            % improvement too small: EXIT
        end
        cost = costnew;
        
    else
        % increase lambda (get closer to gradient descent)
        lambda = lambda * lambdaFactor;
        
        if lambda>lambdaMax && iter > minIter
            break;            % lambda too large: EXIT
        end
    end
end

% print final result if necessary
if flgPrint==2
    fprintf('Iterations = %d;  Cost = %.6f\n', iter, cost);
end

if flgFix
    warning('ilqg_det had to limit large numbers, results may be inaccurate');
end
% close(waithandle);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulate controlled system, compute trajectory and cost
function [x, cost] = simulate(fnDyn, fnCost, dt, x0, u, maxValue)

szX = size(x0,1);
n = size(u,2)+1;

% Initialize simulation
x = zeros(szX, n);
x(:,1) = x0;
cost = 0;

% Run simulation with substeps
for k = 1:n-1
    f = fnDyn(x(:,k), u(:,k));    
    x(:,k+1) = x(:,k) + dt * f;    
    x(:,k+1) = fixBig(x(:,k+1), maxValue);
    cost = cost + dt * fnCost(x(:,k), u(:,k), k);
end

% Adjust for final cost, subsample trajectory
cost = cost + fnCost(x(:,end), NaN, NaN);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute optimal control law
function [l, L] = uOptimal(g, G, H, u, uMin, uMax, lambda)

% Eigenvalue decomposition, modify eigenvalues
if isnan(H)
    H = zeros(size(H));
end
[V,D] = eig(H);

% Ignore imaginary values
if ~isreal(V)
    V = real(V);
end
if ~isreal(D)
    D = real(D);
end

d = diag(D);
d(d<0) = 0;
d = d + lambda;

% Inverse modified Hessian, unconstrained control law
H1 = V*diag(1./d)*V';

l = -H1*g;
L = -H1*G;

% Enforce constraints
l = min(max(l+u, uMin), uMax) - u;

% Modify L to reflect active constraints
L((l+u<=uMin)|(l+u>=uMax), :) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Limit numbers that have become suspiciously large
function s = fixBig(s, maxValue)
global flgFix;

if ~isreal(s)
    s = real(s);
end

ibad = (abs(s)>maxValue);
s(ibad) = sign(s(ibad)) * maxValue;

if any(ibad(:))
    flgFix = 1;
end
