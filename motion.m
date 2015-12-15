function [ out ] = motion( IC ,t )
%Problem 3 Hw 7 
%   Plots movement over time evolution
m1 = 4;
m2 = 3;

theta(1) = IC(1);
l = IC(2);
g = IC(3);
thetad(1) = IC(4);

out(1) = thetadd;
out(2) = -(981*sin(thetad))/10;
out(3) = xdd;
out(4) = -(2943*cos(theta)*sin(theta))/(10*(3*cos(theta)^2 + 70));

end

