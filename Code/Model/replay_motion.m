function replay_motion(x, target, step, sample,flg)

parameters;  % Load arm parameters

for i=1:sample
    th1=x(1,i);th2=x(2,i);
    T00=[cos(th1),-sin(th1),0;sin(th1),cos(th1),0;0,0,1];
    T01=[cos(th2),-sin(th2),l1_;sin(th2),cos(th2),0;0,0,1];
    T12=[1,0,l2_;0,1,0;0,0,1];
    T02=T00*T01*T12;
     
    xy1=[1,0,0;0,1,0]*T00*T01*[0;0;1];
    xy2=[1,0,0;0,1,0]*T02*[0;0;1];
    xy=[[0;0],xy1,xy2];
    
    clf;
    hold on;
    plot(target(1,1),target(2,1),'ro');    
    plot(xy(1,:), xy(2,:), 'ko-','linewidth',2);
    
    if strcmp(flg,'tra') == 1    % Display trajectory
        end_effector(i,:) = xy2;
        plot(end_effector(:,1),end_effector(:,2),'b-','linewidth',1);        
    elseif strcmp(flg,'hist') == 1    % Display post arm postures
        if i == 1
            state1(1,:) = xy(1,:);
            state2(1,:) = xy(2,:);
        end
        interval = 7;
        if rem(i,interval) == 0
            state1(i/interval+1,:) = xy(1,:);
            state2(i/interval+1,:) = xy(2,:);
        end
        for j = 1:size(state1,1)
            plot(state1(j,:), state2(j,:), 'o-','linewidth',1);
        end
    end
    hold off;    
    axis equal;
    axis([-0.8,0.8,-0.2,1.0]);
    grid;
    title(sprintf('%5.2f',(i-1)*step));    
    drawnow;
end
