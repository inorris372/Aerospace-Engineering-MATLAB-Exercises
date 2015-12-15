%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Lab 4 Report Problem      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p = [-3 -2 -1 1 2 3];
dl = [4.625 5.75 6.75 9.375 10.5 11.25];
dr = [9.5 8.25 7 4 2.75 1.75];
theta = size(p);
tau = zeros(1497);
s = size(tau);
a = size(s);
I = 0.52147; %inches^4
t = .08; %inches
y = 1.2505; %inches
V = .042243; %kg
b = size(s);

for x = 1:6
    theta(x) = asin((dl(x)-dr(x))/17.625);
end

for s = 1:size(tau)
    a(s) = (s-1)/1000;
    tau(s) = V/(I*t)*y*t*(s./1000);
    b(s) = (s-1)*2.51/1497;
end



figure(1)
plot(theta,p)
title('Twist Angle (x-axis) vs Location of Weight (y-axis)')

figure(2)
plot(tau,a)
title('Shear Stress Distribution on Flanges (tau = x-axis, s = y-axis)')

figure(3)
plot(b,tau(1497))
title('Shear Stress Distribution on Web (tau = y-axis, s = x-axis)')




