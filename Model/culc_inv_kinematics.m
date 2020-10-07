function out = culc_inv_kinematics(xy)

parameters;  % Load arm parameters

x = xy(1,:);
y = xy(2,:);
l0_ = sqrt(x.^2 + y.^2);

cos_theta2 = -((x.^2+y.^2)-(l1_.^2+l2_.^2))./(2.*l1_.*l2_);
sin_theta2 = sqrt((2.*l1_.*l2_).^2-(l0_.^2 - (l1_.^2+l2_.^2)).^2)./(2.*l1_.*l2_);

theta2 = atan2(real(sin_theta2),real((-cos_theta2)));

kc = l1_ + l2_.*cos(theta2);
ks = l2_.*sin(theta2);

cos_theta1 = (kc.*x+ks.*y)./(kc.^2+ks.^2);
sin_theta1 = (-ks.*x+kc.*y)./(kc.^2+ks.^2);

theta1 = atan2(sin_theta1,cos_theta1);

out = [theta1; theta2];