function [ ysdot ] = EquationsOfMotion(t,s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global mu
ysdot = zeros(4,1);  % a column vector
x(1) = s(1) ;
dxbydt(1) = s(3);
y(1) = s(2);
dybydt(1) = s(4);
% z(0) = 0;
% dzbydt(0) = 0;

y1s = x;
y2s = dxbydt;
y3s = y;
y4s = dybydt;

r1 = sqrt((y1s+mu)^2+y3s^2);
r2 = sqrt((y1s-1+mu)^2+y3s^2);

ysdot(1) = y2s;
ysdot(2) = 2*y4s + y1s -(1-mu)*(y1s+mu)/r1^3-mu*(y1s-1+mu)/r2^3;
ysdot(3) = y4s;
ysdot(4) = -2*y2s + y3s - (1-mu)*y3s/r1^3-mu*y3s/r2^3;

% %the distances
% r1=sqrt((mu+y(1))^2+(y(2))^2+(y(3))^2);
% r2=sqrt((1-mu-y(1))^2+(y(2))^2+(y(3))^2);
% %masses
% m1=1-mu;
% m2=mu;
% 
% ydot=[y(4); 
%     y(5); 
%     y(6); 
%     y(1)+2*y(5)+G*m1*(-mu-y(1))/(r1^3)+G*m2*(1-mu-y(1))/(r2^3); 
%     y(2)-2*y(4)-G*m1*(y(2))/(r1^3)-G*m2*y(2)/(r2^3); 
%     -G*m1*y(3)/(r1^3)-G*m2*y(3)/(r2^3)];
end

