function [xy, xydot, xydotdot, xydotdot_theta1, xydotdot_theta2] = culc_kinematics(theta)

parameters;  % Load arm parameters

theta1 = theta(1, :);
theta2 = theta(2, :);

xy = [l1_*cos(theta1) + l2_*cos(theta1 + theta2); ...
    l1_*sin(theta1) + l2_*sin(theta1 + theta2)];

xydot = [-l1_*sin(theta1)-l2_*sin(theta1 + theta2), -l2_*sin(theta1 + theta2);...
            l1_*cos(theta1)+l2_*cos(theta1 + theta2), l2_*cos(theta1 + theta2)];
        
xydotdot = [-l1_*cos(theta1)-l2_*cos(theta1 + theta2), -l2_*cos(theta1 + theta2);...
            -l1_*sin(theta1)+l2_*sin(theta1 + theta2), -l2_*sin(theta1 + theta2)];
        
xydotdot_theta1 = [-l1_*cos(theta1)-l2_*cos(theta1 + theta2), -l2_*cos(theta1 + theta2); ...
            -l1_*sin(theta1)-l2_*sin(theta1 + theta2), -l2_*sin(theta1 + theta2)];

xydotdot_theta2 = [-l2_*cos(theta1 + theta2), -l2_*cos(theta1 + theta2);...
            -l2_*sin(theta1 + theta2), -l2_*sin(theta1 + theta2)];
        

