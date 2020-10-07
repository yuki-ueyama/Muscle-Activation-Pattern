%% Link parameters
% Mass [kg]
m1_   = 0.3;
m2_   = 0.3;

% Length [m]
l1_ = 0.15;
l2_ = 0.21;

% Center of the mass [m]
s1_ = 0.07; % It does not affect
s2_ = 0.12;

% Moment of inertia [kg/m^2]
i1_ = 5.0*1e-3;
i2_ = 9.0*1e-3;

% Joint viscosity [NÅEmÅEs/rad]
b11_  = 5.0*1e-3;
b22_  = 5.0*1e-3;
b12_  = 2.5*1e-3;
b21_  = 2.5*1e-3;

%% Muscle parameters
% % Moment arm [m]
% J = [2.0 -2.0 0 0 1.5 -2.0;
%     0 0 2.0 -2.0 2.0 -1.5]./100;
% 
% % Psiological cross-sectional area [cm^2]
% PCSA = [11; 7.0; 7.0; 10.5; 7.0; 12]; 
% 
% % Optimal muscle length [m]
%  L0 = [8.0; 5.0; 10.0; 4.5; 6.5; 4.0]./100;
 
% Moment arm [m]
J = [1.5 -1.5 0 0 1.5 -1.5;
    0 0 1.5 -1.5 1.5 -1.5]./100;

% Psiological cross-sectional area [cm^2]
PCSA = [10; 10; 10; 10; 10; 10];

% Optimal muscle length [m]
L0 = [8; 8; 8; 8; 8; 8]./100;

% Absolute muscle force [N/cm^2]
Fm = 32;
 
% Optimal joint angle [rad]
eq_theta = [15 5 0 0 15 5;
    0 0 90 110 100 100]*pi/180;

% Time constants of muscle dynamics [s]
t_act = 0.050;
t_deact = 0.066;