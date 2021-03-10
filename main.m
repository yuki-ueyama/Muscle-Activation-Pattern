function out = main(varargin)

addpath('Model');   % Add sub-folder

%% Simulation setting ----------------------------------------------
% simset.ctl: control mode
%   t: terminal control;
%   s: stabilization control; 
%   c: constraint states at a temporal time;
% simset.cost: structure of the cost function
%   p: postion;
%   pv: postion + velocity;
%   pvf: postion + velocity + force(or torque)
% simset.coor: coordinate of the cost funtion
%   c: Cartesian (endpoint);
%   j: joint
simset.ctl  = 't'; 
simset.cost = 'pvf'; 
simset.coor = 'c';

% Scaling factor for sensitive analysis
if nargin == 1
    simset.scale = varargin{1};
elseif nargin < 1
     simset.scale = 1;
else
    disp('Invalid inputs of main !!!');
end

% Set parameters -------------------------------------------
if strcmp(simset.ctl, 't')
    simset.mt = 0.40;      % Movement time [s]
    filename_opt = '';
elseif strcmp(simset.ctl, 's')
    simset.mt = 0.50;
    filename_opt = 'stab_';
elseif strcmp(simset.ctl, 'c')
    simset.mt = 0.50;
    filename_opt = 'tmp_';
else
    disp('Error; Set simset.ctl to ''s'', ''t'' or ''c'' !!!');
    return;
end

simset.dt = 0.005;       % Size of time-step [s]
N = simset.mt/simset.dt + 1;    % Number of time steps
time = simset.dt*(0:N-1);

parameters;  % Load parameters

distance = 0.08;           % Reaching distance [m]
directions = 32;            % Number of reaching directions

% Initial hand position
xy0 = [0.08; 0.22];      % Almost equal to a muscle equilibrium point

% Target positions
alpha = 2*pi*(0:directions-1)./directions;  % Reaching directions
hand_target = repmat(xy0, [1 directions]) + distance*[cos(alpha); sin(alpha)];
% joint_target = culc_inv_kinematics(hand_target);

%  Initialize state and control signal
initial = culc_inv_kinematics(xy0);
x0 = [initial; zeros(8,1)];
uMin = 0.0; uMax = 1.0; % Limitation of control signal
u0 = 0.1*ones(6, 1);    % Default

% Variables -------------------------------------------------------
x = zeros(10, N, directions);           % State vector
u = zeros(6, N-1, directions);          % Control signal
torque = zeros(2, N, directions);       % Joint torque
hand_pos = zeros(2, N, directions);     % Hand position
hand_vel = zeros(2, N, directions);     % Hand velocity
hand_force = zeros(2, N, directions);   % Hand force

%% Loop start ------------------------------------------------------
for d = 1:directions
    
    % Cost function
     fnCost = @(x_,u_,t_) cost_function(x_,u_,t_, hand_target(:, d), simset);
   
    % Arm dynamics function
    fnDyn = @(x_,u_) dynamics_function(x_, u_);
    
    % Compute optimal control using Iterative LQG
    [x(:,:,d), u(:,:,d), cost] = ...
        ilqg_det(fnDyn, fnCost, simset.dt, N, x0, u0, uMin, uMax);

    torque(:, :, d) = J*muscle_tension(x(:, :, d)); % Joint torque
    
    % Caluculate forward kinematics
    hand_pos(:, :, d) = culc_kinematics(x(1:2, :, d)) - repmat(xy0, [1 N]);    % Hand position
    Jacob = squeeze(culc_jacobian(x(1:2, :, d)));                              % Jacobian matrix
    for k=1:N
        hand_vel(:, k, d) = squeeze(Jacob(:, :, k))*x(3:4, k, d); % Hand velocity
        hand_force(:, k, d) = squeeze(Jacob(:, :, k))'\squeeze(torque(:, k, d));    % Hand force
    end
    
    %% Plot functions ----------------------------------------------
    % Plot result
    figure(1);  
    subplot(2,2,1); % Hand trajectory
    plot(0, 0, 'ro');
    hold on;
    plot(hand_target(1, 1:d)-repmat(xy0(1, 1), [1 d]), ...
        hand_target(2, 1:d)-repmat(xy0(2, 1), [1 d]),'go');
    plot(squeeze(hand_pos(1, :, 1:d)), squeeze(hand_pos(2, :, 1:d)));
    xlim ([-0.1 0.1]);
    ylim ([-0.1 0.1]);
    xlabel('x-position [m]');
    ylabel('y-position [m]');
    axis square;
    hold off;
    subplot(2,2,2) % Hand velocity
    plot(time, squeeze(sqrt(sum(hand_vel(1:2, :, 1:d).^2)))');
    xlabel('Time [s]');
    ylabel('Speed [m/s]');
    axis square;
    subplot(2,2,3) % Hand force
    plot(time, squeeze(sqrt(sum(hand_force(1:2, :, 1:d).^2)))');
    xlabel('Time [s]');
    ylabel('Force norm [N]');
    axis square;
    subplot(2,2,4) % Joint torque
    plot(time, squeeze(torque(1, :, 1:d))', '-'); hold on;
    plot(time, squeeze(torque(2, :, 1:d))', ':');
    xlabel('Time [s]');
    ylabel('Torque [Nm]');
    axis square;
    hold off;
    
    % Replay limb motion
    figure(2);
    replay_motion(x(:, : , d), hand_target(:, d), simset.dt, N, 'hist');
    
    drawnow;
end
%% Loop end --------------------------------------------------------

%% Plot muscle activation patterns ---------------------------------
figure(3)
colormap('default');    % colormap('gray');
angles = ((0:directions)*360/directions);
U = zeros(6, N, directions+1);
U(:,:,1:directions) = x(5:10,:,:);
U(:,:,directions+1) = U(:,:,1);
subplot(2,3,1);
imagesc(time, angles, squeeze(U(1,:,:))');
title('SF');
ylim([0 360]);
set(gca,'ytick',[0 180 360]);
subplot(2,3,2);
imagesc(time, angles, squeeze(U(5,:,:))');
title('BF');
ylim([0 360]);
set(gca,'ytick',[0 180 360]);
subplot(2,3,3);
imagesc(time, angles, squeeze(U(3,:,:))');
title('EF');
ylim([0 360]);
set(gca,'ytick',[0 180 360]);
subplot(2,3,4);
imagesc(time, angles, squeeze(U(2,:,:))');
title('SX');
ylabel('Target direction [deg]');
xlabel('Time [s]');
ylim([0 360]);
set(gca,'ytick',[0 180 360]);
subplot(2,3,5);
imagesc(time, angles, squeeze(U(6,:,:))');
title('BX');
ylim([0 360]);
set(gca,'ytick',[0 180 360])
subplot(2,3,6);
imagesc(time, angles, squeeze(U(4,:,:))');
title('EX');
ylim([0 360]);
set(gca,'ytick',[0 180 360]);
%     colorbar('Ticks',[0,0.5,1])

%% Save result -----------------------------------------------------  
filename = ['data_' filename_opt simset.cost int2str(directions)];
if nargin == 1
    filename = [filename 'x' num2str(simset.scale) '.mat'];
else
    filename = [filename '.mat'];
end
save(filename, 'x', 'u', 'x0', 'xy0', 'torque', 'hand_force', 'hand_pos', ...
    'hand_target', 'hand_vel', 'N', 'directions', 'time', 'simset');
out = filename;

