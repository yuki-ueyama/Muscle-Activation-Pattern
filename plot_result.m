clear;
clf;

dname = 'Data';   % Directory name
% % File names of data
% fname = {'data_p32.mat', 'data_pv32.mat', 'data_pf32.mat', 'data_pvf32.mat'};
fname = {'data_tmp_p32.mat', 'data_tmp_pv32.mat', 'data_tmp_pf32.mat', 'data_tmp_pvf32.mat'};
fname2 = {'data_stab_p32.mat', 'data_stab_pv32.mat', 'data_stab_pf32.mat', 'data_stab_pvf32.mat'};

% dname = 'SensitivityAnalysis';   % Directory name
% % File names of data
% fname = {'data_p32x0.001.mat', 'data_pv32x0.001.mat', 'data_pf32x0.001.mat', 'data_pvf32x0.001.mat'};
% fname2 = {'data_stab_p32x0.001.mat', 'data_stab_pv32x0.001.mat', 'data_stab_pf32x0.001.mat'};
% fname = {'data_p32x0.01.mat', 'data_pv32x0.01.mat', 'data_pf32x0.01.mat', 'data_pvf32x0.01.mat'};
% fname2 = {'data_stab_p32x0.01.mat', 'data_stab_pv32x0.01.mat', 'data_stab_pf32x0.01.mat'};
% fname = {'data_p32x0.1.mat', 'data_pv32x0.1.mat', 'data_pf32x0.1.mat', 'data_pvf32x0.1.mat'};
% fname2 = {'data_stab_p32x0.1.mat', 'data_stab_pv32x0.1.mat', 'data_stab_pf32x0.1.mat'};
% fname = {'data_p32x10.mat', 'data_pv32x10.mat', 'data_pf32x10.mat', 'data_pvf32x10.mat'};
% fname2 = {'data_stab_p32x10.mat', 'data_stab_pv32x10.mat', 'data_stab_pf32x10.mat'};
% fname = {'data_p32x100.mat', 'data_pv32x100.mat', 'data_pf32x100.mat', 'data_pvf32x100.mat'};
% fname2 = {'data_stab_p32x100.mat', 'data_stab_pv32x100.mat', 'data_stab_pf32x100.mat'};
% fname = {'data_p32x1000.mat', 'data_pv32x1000.mat', 'data_pf32x1000.mat', 'data_pvf32x1000.mat'};
% fname2 = {'data_stab_p32x1000.mat', 'data_stab_pv32x1000.mat', 'data_stab_pf32x1000.mat'};

% dname = 'SensitivityAnalysis';   % Directory name
% File names of data
% fname = {'data_tmp_p32x0.001.mat', 'data_tmp_pv32x0.001.mat', 'data_tmp_pf32x0.001.mat', 'data_tmp_pvf32x0.001.mat'};
% fname2 = {'data_stab_p32x0.001.mat', 'data_stab_pv32x0.001.mat', 'data_stab_pf32x0.001.mat'};
% fname = {'data_tmp_p32x0.01.mat', 'data_tmp_pv32x0.01.mat', 'data_tmp_pf32x0.01.mat', 'data_tmp_pvf32x0.01.mat'};
% fname2 = {'data_stab_p32x0.01.mat', 'data_stab_pv32x0.01.mat', 'data_stab_pf32x0.01.mat'};
% fname = {'data_tmp_p32x0.1.mat', 'data_tmp_pv32x0.1.mat', 'data_tmp_pf32x0.1.mat', 'data_tmp_pvf32x0.1.mat'};
% fname2 = {'data_stab_p32x0.1.mat', 'data_stab_pv32x0.1.mat', 'data_stab_pf32x0.1.mat'};
% fname = {'data_tmp_p32x10.mat', 'data_tmp_pv32x10.mat', 'data_tmp_pf32x10.mat', 'data_tmp_pvf32x10.mat'};
% fname2 = {'data_stab_p32x10.mat', 'data_stab_pv32x10.mat', 'data_stab_pf32x10.mat'};
% fname = {'data_tmp_p32x100.mat', 'data_tmp_pv32x100.mat', 'data_tmp_pf32x100.mat', 'data_tmp_pvf32x100.mat'};
% fname2 = {'data_stab_p32x100.mat', 'data_stab_pv32x100.mat', 'data_stab_pf32x100.mat'};
% fname = {'data_tmp_p32x1000.mat', 'data_tmp_pv32x1000.mat', 'data_tmp_pf32x1000.mat', 'data_tmp_pvf32x1000.mat'};
% fname2 = {'data_stab_p32x1000.mat', 'data_stab_pv32x1000.mat', 'data_stab_pf32x1000.mat'};
%% Plot simulated results
% figure(1);
figure('Name','Trajectories', 'Position', [0 100 800 800]);
nset = length(fname);
umax = zeros(nset, 6);
amax = zeros(nset, 6);

for n=1:nset
    load( [dname '\' fname{n}]);    
    
    if directions == 32
        d_list = 1:2:directions;
    else
        d_list = 1:directions;
    end
    
    % Hand trajectory
    subplot(nset,4,4*(n-1)+1);
   h = plot(100*hand_target(1, d_list)-100*repmat(xy0(1, 1), [1 length(d_list)]), ...
        100*hand_target(2, d_list)-100*repmat(xy0(2, 1), [1 length(d_list)]),'o', 'Color', [0.8 0.8 0.8]);
    set(h, 'MarkerFaceColor', get(h,'Color')); 
    hold on;
    plot(100*squeeze(hand_pos(1, :, d_list)), 100*squeeze(hand_pos(2, :, d_list)));
    xlim ([-10 10]);
    ylim ([-10 10]);
    if n == 1
        title('Hand path');
    end
    xlabel('x-position [cm]');
    ylabel('y-position [cm]');
    axis square; box off; hold off;
    
    % Hand velocity
    subplot(nset,4,4*(n-1)+2)
    plot(time, squeeze(sqrt(sum(hand_vel(1:2, :, d_list).^2)))');
    ylim ([0 simset.mt]);
    if n == 1
        title('Hand speed');
    end
    xlabel('Time [s]');
    ylabel('Velocity [m/s]');
    xticks([0 0.2 0.4 0.5]);
    axis square; box off;
    
    % Hand force
    subplot(nset,4,4*(n-1)+3)
    plot(time, squeeze(sqrt(sum(hand_force(1:2, :, d_list).^2)))');
     ylim ([0 4.0]);
%     ylim ([0 3.0]);
    if n == 1
        title('Hand force');
    end
    xlabel('Time [s]');
    ylabel('Force norm [N]');
    xticks([0 0.2 0.4 0.5]);
    axis square; box off;
    
    % Joint torque
    subplot(nset,4,4*(n-1)+4)
    plot(time, squeeze(torque(1, :, d_list))', '-'); hold on;
    plot(time, squeeze(torque(2, :, d_list))', ':');
    ylim ([-0.5, 0.5]);
    if n == 1
        title('Joint torque');
    end
    xlabel('Time [s]');
    ylabel('Torque [Nm]');
    xticks([0 0.2 0.4 0.5]);
    axis square; box off; hold off;
    
    for k=1:6
        umax(n,k) = max(max(u(k,:,:)));
        amax(n,k) = max(max(x(k+4,:,:)));
    end
end

%% Plot muscle activation patterns
% figure('Name','Muscle activation patterns', 'Position', [0 100 400 400]);
colormap('default');    % colormap('gray');
thres = 0;
clims = [repmat(thres, [6 1]) max(amax)'];
angles = ((0:directions)*360/directions);
U = zeros(6, N, directions+1, nset);

for n=1:nset
    load( [dname '\' fname{n}]);
        
    if directions == 32
        d_list = 1:2:directions;
    else
        d_list = 1:directions;
    end
        
    %     U(:,2:N,1:directions, n) = u;
    U(:,1:N,1:directions, n) = x(5:10,:,:);
    U(:,:,directions+1, n) = U(:,:,1, n);
    f = figure(n+1);
    set(f,'Position', [100 100 350 250]);
    subplot(2,3,1);
    imagesc(time, angles, squeeze(U(1,:,:, n))', clims(1,:));
    title('SF');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    subplot(2,3,2);
    imagesc(time, angles, squeeze(U(5,:,:, n))', clims(5,:));
    title('BF');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    subplot(2,3,3);
    imagesc(time, angles, squeeze(U(3,:,:, n))', clims(3,:));
    title('EF');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    subplot(2,3,4);
    imagesc(time, angles, squeeze(U(2,:,:, n))', clims(2,:));
    title('SX');
    ylabel('Target direction [deg]');
    xlabel('Time [s]');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    subplot(2,3,5);
    imagesc(time, angles, squeeze(U(6,:,:, n))', clims(6,:));
    title('BX');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360])
    subplot(2,3,6);
    imagesc(time, angles, squeeze(U(4,:,:, n))', clims(4,:));
    title('EX');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    %     colorbar('Ticks',[0,0.5,1])
end

% % Plot exapmle muscle activations in a direction
% figure(100)
% tmp = 150;
% subplot(2,3,1)
% plot(time, squeeze(U(1,:,round(tmp/360*directions), :))');
% box off;
% ylim([0 0.4]);
% subplot(2,3,2)
% plot(time, squeeze(U(5,:,round(tmp/360*directions), :))');
% box off;
% ylim([0 0.4]);
% subplot(2,3,3)
% plot(time, squeeze(U(3,:,round(tmp/360*directions), :))');
% box off;
% ylim([0 0.4]);
% subplot(2,3,4)
% plot(time, squeeze(U(2,:,round(tmp/360*directions), :))');
% box off;
% ylim([0 0.4]);
% subplot(2,3,5)
% plot(time, squeeze(U(6,:,round(tmp/360*directions), :))');
% box off;
% ylim([0 0.4]);
% subplot(2,3,6)
% plot(time, squeeze(U(4,:,round(tmp/360*directions), :))');
% box off;
% ylim([0 0.4]);

% Plot extra case
nset = length(fname2);
umax = zeros(nset, 6);
amax = zeros(nset, 6);

for n=1:nset
    load( [dname '\' fname2{n}]);
    U(:,1:N,1:directions, n) = x(5:10,:,:);
    U(:,:,directions+1, n) = U(:,:,1, n);
    % Hand trajectory
    f = figure(10+n);    
    set(f,'Position', [100 100 400 400]);
    subplot(2,2,1);
    h = plot(100*hand_target(1, d_list)-100*repmat(xy0(1, 1), [1 length(d_list)]), ...
        100*hand_target(2, d_list)-100*repmat(xy0(2, 1), [1 length(d_list)]),'o', 'Color', [0.8 0.8 0.8]);
    set(h, 'MarkerFaceColor', get(h,'Color')); 
    hold on;
    plot(100*squeeze(hand_pos(1, :, d_list)), 100*squeeze(hand_pos(2, :, d_list)));
    xlim ([-10 10]);
    ylim ([-10 10]);
    title('Hand path');
    xlabel('x-position [cm]');
    ylabel('y-position [cm]');
    axis square; box off; hold off;
    % Hand velocity
    subplot(2,2,2)
    plot(time, squeeze(sqrt(sum(hand_vel(1:2, :, d_list).^2)))');
    ylim ([0 0.4]);
    xlim ([0 0.5]);
    title('Hand speed');
    xlabel('Time [s]');
    ylabel('Velocity [m/s]');
    xticks([0 0.2 0.4 0.5]);
    axis square; box off;
    % Hand force
    subplot(2,2,3)
    plot(time, squeeze(sqrt(sum(hand_force(1:2, :, d_list).^2)))');
    ylim ([0 4]);
    xlim ([0 0.5]);
    xticks([0 0.2 0.4 0.5]);
    title('Hand force');
    xlabel('Time [s]');
    ylabel('Force norm [N]');
    axis square; box off;
    % Joint torque
    subplot(2,2,4)
    plot(time, squeeze(torque(1, :, d_list))', '-'); hold on;
    plot(time, squeeze(torque(2, :, d_list))', ':');
    xlim ([0 0.5]);
    ylim ([-0.5 0.5]);
    xticks([0 0.2 0.4 0.5]);
    title('Joint torque');
    xlabel('Time [s]');
    ylabel('Torque [Nm]');
    axis square; box off; hold off;
    
    f = figure(20+n);
    set(f,'Position', [100 100 350 250]);
    subplot(2,3,1);
    imagesc(time, angles, squeeze(U(1,:,:, n))');
    title('SF');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    subplot(2,3,2);
    imagesc(time, angles, squeeze(U(5,:,:, n))');
    title('BF');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    subplot(2,3,3);
    imagesc(time, angles, squeeze(U(3,:,:, n))');
    title('EF');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    subplot(2,3,4);
    imagesc(time, angles, squeeze(U(2,:,:, n))');
    title('SX');
    ylabel('Target direction [deg]');
    xlabel('Time [s]');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
    subplot(2,3,5);
    imagesc(time, angles, squeeze(U(6,:,:, n))');
    title('BX');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360])
    subplot(2,3,6);
    imagesc(time, angles, squeeze(U(4,:,:, n))');
    title('EX');
    ylim([0 360]);
    set(gca,'ytick',[0 180 360]);
end