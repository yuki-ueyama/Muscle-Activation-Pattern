function [T, Tdot_a, Tdotdot_a, Tdot_theta, Tdotdot_theta, ...
    Tdot_thetadot, Tdotdot_thetadot] = muscle_tension(x)

szX2 = size(x, 2);
szM = 6; % Number of muscles

parameters;  % Load arm parameters

% Muscle length
L = zeros(szM, szX2);
for k = 1:szX2
    L(:, k) = diag(J'*(eq_theta - repmat(x(1:2, k), [1 szM])));
end

% Muscle shortening velocity
L_d = -J' * x(3:4, :) ;

% Normalize
l = ones(szM, szX2) + L./repmat(L0, [1 szX2]); 
l_v = L_d./(repmat(L0, [1 szX2])); 

% Compute muscle tensions
[T, Tdot_a, Tdotdot_a, Tdot_l, Tdotdot_l, Tdot_ldot, Tdotdot_ldot] ...
    = muscle_model(x(5:10, :), l, l_v, PCSA, Fm);

if nargout>1
    Tdot_theta = diag(L0.*Tdot_l)*J';
    Tdot_thetadot = diag(L0.*Tdot_ldot)*J';
    
    Tdotdot_theta = zeros(2, 2, szM);
    Tdotdot_thetadot = zeros(2, 2, szM);
    
    tmp = diag(L0.*Tdotdot_l)*J';
    tmp2 = diag(L0.*Tdotdot_ldot)*J';
    
    for k = 1:szM
        Tdotdot_theta(:, :, k) = diag(tmp(k,:));
        Tdotdot_thetadot(:, :, k) = diag(tmp2(k,:));
    end    
end
