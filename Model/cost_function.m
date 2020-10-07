function [l, l_x, l_xx, l_u, l_uu, l_ux] = cost_function(x, u, t, target, simset)
if isnan(u)
    u = 0;
end

szX = size(x,1);
szU = size(u,1);
szX2 = size(x,2);
target = repmat(target,[1,szX2]);

% Set cost weights
wp = simset.scale*1.0*1e+3;  % Positional error
wv = simset.scale*1.0*1e+2;  % Terminal velocity
wf = simset.scale*1.0*1e+1;  % Ternimal force
wu = 1.0;                    % Effort

if strcmp(simset.cost, 'p')
    wv = 0; wf = 0;
elseif strcmp(simset.cost, 'pv')
    wf = 0;
elseif strcmp(simset.cost, 'pf')
    wv = 0;
elseif strcmp(simset.cost, 'pvf') 
else
    disp('Error; Set simset.cost adequately !!!');
end

if strcmp(simset.ctl, 't')
    N = simset.mt/simset.dt + 1;  % final step
elseif strcmp(simset.ctl, 's')
    stab_period = 0.1;                                      % Stabilization period: 0.1 [s]
    N = simset.mt/simset.dt + 1 - stab_period/simset.dt;    % Stabilization step
    wp = wp/stab_period;
    wv = wv/stab_period;
    wf = wf/stab_period;
else
    disp('Error; Set simset.ctl to ''s'' or ''t'' !!!');
end

l = zeros(1,szX2);
vel = zeros(2,szX2);
force = zeros(2,szX2);

% Compute cost
if isnan(t) || t > N
    [e, edot, edotdot] = culc_kinematics(x(1:2,:));
    [Jacob, Jdot, Jdot_theta1, Jdot_theta2, Jdotdot_theta1, Jdotdot_theta2] = culc_jacobian(x(1:2,:));
    [T, Tdot_a, Tdot_aa, Tdot_theta, Tdotdot_theta, ...
        Tdot_thetadot, Tdotdot_thetadot] = muscle_tension(x);
    
    parameters;  % Load parameters
    torque = J*T;
    
    for k=1:szX2
        vel(:, k) = squeeze(Jacob(:,:,k))*x(3:4, k);
        force(:, k) = squeeze(Jacob(:,:,k))'\torque(:, k);
    end
             
    if strcmp(simset.coor, 'c')
        l(1,:) = wp*sum((e-target).^2) + wv*sum(vel.^2) + wf*sum(force.^2);
    elseif strcmp(simset.coor, 'j')
        joint_target = culc_inv_kinematics(target);
        l(1,:) = wp*sum((x(1:2,:)-joint_target).^2) + wv*sum(x(3:4,:).^2) + wf*sum(torque.^2);
    else
        disp('Error; Set mode.coor to ''c'' or ''j'' !!!');
    end  

    l(1,:) = l(1,:) + wu*sum(u.^2);
    
else % Effort cost
    l(1,:) =  wu*sum(u.^2) ;
end

% Compute derivatives of cost
if nargout > 1
    l_u =  2*wu*u;
    l_uu = 2*wu*eye(szU);
    l_ux = zeros(szU, szX);
    l_x = zeros(szX, 1);
    l_xx = zeros(szX, szX);
    l_e = zeros(6, 1);
    e_x = zeros(6, szX);
    
    if isnan(t) || t >N
        if strcmp(simset.coor, 'c')           % Cartesian coordinate
            l_e(1:2, 1) = 2*wp*(e(:,1)-target(:,1)) ;
            l_e(3:4, 1) = 2*wv*vel(:,1) ;
            l_e(5:6, 1) = 2*wf*force(:,1);
            
            e_x(1:2, 1:2) = edot;
            e_x(3:4, 1) = Jdot_theta1*x(3:4,:);
            e_x(3:4, 2) = Jdot_theta2*x(3:4,:);
            e_x(3:4, 3:4) = Jacob;
            e_x(5:6, 1:2) = Jacob'\(J*Tdot_theta);
            e_x(5:6, 1) = e_x(5:6, 1) - (Jacob'\Jdot_theta1'/Jacob')*torque(:,1);
            e_x(5:6, 2) = e_x(5:6, 2) - (Jacob'\Jdot_theta2'/Jacob')*torque(:,1);
            e_x(5:6, 3:4) = Jacob'\(J*Tdot_thetadot);
            e_x(5:6, 5:10) = Jacob'\(J*diag(Tdot_a));
                        
            l_x = e_x'*l_e;
            
            l_ee(1:2) = 2*wp*ones(1,2);
            l_ee(3:4) = 2*wv*ones(1,2);
            l_ee(5:6) = 2*wf*ones(1,2);
            
            tmp = zeros(2,2);
            tmp2 = zeros(2,2);
            for k = 1:6
                tmp = tmp + diag(J(:, k))*squeeze(Tdotdot_theta(:,:,k));
                tmp2 = tmp2 +diag(J(:, k))*squeeze(Tdotdot_thetadot(:,:,k));
            end
            
            e_xx(1:2, 1:2) = edotdot;
            e_xx(3:4, 1:2) = Jdot;
            e_xx(3:4, 1) =  e_xx(3:4, 1) + Jdotdot_theta1*x(3:4,:);
            e_xx(3:4, 2) =  e_xx(3:4, 2) + Jdotdot_theta2*x(3:4,:);
            e_xx(5:6, 1:2) = Jacob'\tmp ...
                - (Jacob'\Jdot'/Jacob')\(J*Tdot_theta) ; 
            e_xx(5:6, 1) = e_xx(5:6, 1) ...
                +  (Jacob'\Jdot_theta1'/Jacob')*(Jdot_theta1'/Jacob')*torque(:,1) ...
                - (Jacob'\Jdotdot_theta1'/Jacob')*torque(:,1) ...
                + (Jacob'\Jdot_theta1')*(Jacob'\Jdot_theta1'/Jacob')*torque(:,1) ...
                +  (Jacob'\Jdot_theta1'/Jacob')*(Jacob'\Jdot_theta1'/Jacob')*torque(:,1); 
            e_xx(5:6, 2) = e_xx(5:6, 2) ...
                +  (Jacob'\Jdot_theta2'/Jacob')*(Jdot_theta2'/Jacob')*torque(:,1) ...
                - (Jacob'\Jdotdot_theta2'/Jacob')*torque(:,1) ...
                + (Jacob'\Jdot_theta2')*(Jacob'\Jdot_theta2'/Jacob')*torque(:,1) ...                
                +  (Jacob'\Jdot_theta2'/Jacob')*(Jacob'\Jdot_theta2'/Jacob')*torque(:,1);   
             e_xx(5:6, 3:4) = Jacob'\tmp2;
            e_xx(5:6, 5:10) = Jacob'\(J*diag(Tdot_aa));
               
            l_xx = e_x'*diag(l_ee)*e_x + diag(e_xx'*l_e);
            
        else            % Joint coordinate
            l_x(1:2) = 2*wp*(x(1:2)-joint_target) + 2*wf*(J*Tdot_theta)'*torque;
            l_x(3:4) = 2*wv*x(3:4) + 2*wf*(J*Tdot_thetadot)'*torque;
            l_x(5:10) = 2*wf*(J*diag(Tdot_a))'*torque;
            
            tmp = zeros(2,2);
            tmp2 = zeros(2,2);
            for k = 1:6
                tmp = tmp + diag(J(:, k))*squeeze(Tdotdot_theta(:,:,k));
                tmp2 = tmp2 +diag(J(:, k))*squeeze(Tdotdot_thetadot(:,:,k));
            end
            l_xx(1:2, 1:2) = 2*wp*eye(2) + 2*wf*((J*Tdot_theta)'^2 ...
                + tmp)*repmat(torque, [1 2]);
            l_xx(3:4, 3:4) = 2*wv*eye(2) + 2*wf*((J*Tdot_thetadot)'^2 ...
                + tmp2)*repmat(torque, [1 2]);
            l_xx(5:10, 5:10) = 2*wf*diag(((J*diag(Tdot_a))'.^2 ...
                + (J*diag(Tdot_aa))')*torque);
        end
    end
end