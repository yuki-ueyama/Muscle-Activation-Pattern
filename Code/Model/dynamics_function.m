function [xdot, xdot_x, xdot_u] = dynamics_function(x, u)

EPS = [1e-7*ones(2, 1); 1e-7*ones(2, 1);1e-7*ones(6, 1)];

parameters;  % Load arm parameters

szX = size(x,1);
szX2 = size(x,2);
szU = size(u,1);

%------------------------ compute inertia I and extra torque H --------
% temp vars
mls = m2_*l1_*s2_;
iml = i1_ + i2_ + m2_*l1_^2;
dd = i2_*iml - i2_^2;

% determinant
det = dd - mls^2*cos(x(2,:)).^2;

% inverse inertia I1
I1_11 = i2_./det;
I1_12 = (-i2_-mls*cos(x(2,:)))./det;
I1_22 = (iml+2*mls*cos(x(2,:)))./det;

% extra torque H (Coriolis, centripetal, friction)
H = [-mls*(2*x(3,:)+x(4,:)).*x(4,:).*sin(x(2,:)) + b11_*x(3,:) + b12_*x(4,:);...
    mls*x(3,:).^2.*sin(x(2,:)) + b22_*x(4,:) + b12_*x(3,:)];

% % ---------------------------------
torque = J*muscle_tension(x);
tau = torque - H ;

theta_dd = [I1_11.*tau(1,:) + I1_12.*tau(2,:); ...
    I1_12.*tau(1,:) + I1_22.*tau(2,:)];

% Muscle activation ------------------------------
t_const = t_deact*ones(szU, szX2);
tmp = u-x(5:10,:) > 0;
t_const(tmp) = t_deact + u(tmp).*(t_act - t_deact);
adot = (u - x(5:10,:))./t_const;

xdot = [x(3:4,:);
    theta_dd;
    adot];
   
%----------- Compute xdot_x using finite differences ------------ 
if nargout>1      
    tmp = diag(ones(1,6)./t_deact);
    tmp2 = (t_deact + u.*(t_act - t_deact)-...
        (u - x(5:10,:)).*(t_act - t_deact))./((t_deact + u.*(t_act-t_deact)).^2);
    tmp3 = u-x(5:10,:)>0;
    tmp(tmp3, tmp3) = diag(tmp2(tmp3)); 
    xdot_u = [zeros(4, szU); tmp];   
    
    x1 = repmat(x, [1,szX]) + diag(EPS);
    x2 = repmat(x, [1,szX]) - diag(EPS);
    uu = repmat(u,[1,szX]);
    
    f1 = dynamics_function(x1, uu);
    f2 = dynamics_function(x2, uu);
    xdot_x = (f1-f2)./2./repmat(EPS, [1 szX]);
else
    xdot_x = [];
    xdot_u =[];
end

