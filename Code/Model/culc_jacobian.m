function [J, Jdot, Jdot_theta1, Jdot_theta2, Jdotdot_theta1, Jdotdot_theta2] = culc_jacobian(x)

parameters; % Load arm parameters
szX3 = size(x,2);
J = zeros(2, 2, szX3);
theta1 = x(1,:);
theta2 = x(2,:);

J(1,1,:) = -l1_*sin(theta1)-l2_*sin(theta1+theta2);
J(1,2,:) = -l2_*sin(theta1+theta2);
J(2,1,:) = l1_*cos(theta1)+l2_*cos(theta1+theta2);
J(2,2,:) = l2_*cos(theta1+theta2);

Jdot(1,1,:) = -l1_*cos(theta1)-l2_*cos(theta1+theta2);
Jdot(1,2,:) = -l2_*cos(theta1+theta2);
Jdot(2,1,:) = -l1_*sin(theta1)-l2_*sin(theta1+theta2);
Jdot(2,2,:) = -l2_*sin(theta1+theta2);

Jdot_theta1(1,1,:) = -l1_*cos(theta1)-l2_*cos(theta1+theta2);
Jdot_theta1(1,2,:) = -l2_*cos(theta1+theta2);
Jdot_theta1(2,1,:) = -l1_*sin(theta1)-l2_*sin(theta1+theta2);
Jdot_theta1(2,2,:) = -l2_*sin(theta1+theta2);

Jdot_theta2(1,1,:) = -l2_*cos(theta1+theta2);
Jdot_theta2(1,2,:) = -l2_*cos(theta1+theta2);
Jdot_theta2(2,1,:) = -l2_*sin(theta1+theta2);
Jdot_theta2(2,2,:) = -l2_*sin(theta1+theta2);

Jdotdot_theta1(1,1,:) = l1_*sin(theta1)+l2_*sin(theta1+theta2);
Jdotdot_theta1(1,2,:) = l2_*sin(theta1+theta2);
Jdotdot_theta1(2,1,:) = -l1_*cos(theta1)-l2_*cos(theta1+theta2);
Jdotdot_theta1(2,2,:) = -l2_*cos(theta1+theta2);

Jdotdot_theta2(1,1,:) = l2_*sin(theta1+theta2);
Jdotdot_theta2(1,2,:) = l2_*sin(theta1+theta2);
Jdotdot_theta2(2,1,:) = -l2_*cos(theta1+theta2);
Jdotdot_theta2(2,2,:) = -l2_*cos(theta1+theta2);
