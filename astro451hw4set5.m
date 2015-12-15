%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ian Norris             %
%     Astrodynamics Hw Set 5  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1

% A sea-level tracking station at 29? north geodetic latitude observes an 
% Earth satellite in transit across its skies with the following data,
%              Local
%  Time (sec)  Sidereal Time (?) Right Ascension (?) Declination (?)
%  0.0         0.0               0.0                 51.5110
%  60.0        0.250684          65.9279             27.9911
%  120.0       0.501369          79.850              14.6609

% Use the Gauss method to find both r and v at the 1-minute point. Ignore 
% iterative improvement. [Hint: |v| ? 8 km/s]
clear;
clc;
phi = 29*pi/180; %Geodetic Latitude
Re = 6378.000000; %Radius of Earth in km
flat = 1/298.6; %Flattening Factor
mu = 3.986e5;
t1 = 0.0;
t2 = 60.0;
t3 = 120.0;
H = 0;

ra1 = 0.0*pi/180;
ra2 = 65.9279*pi/180;
ra3 = 79.850*pi/180;

dec1 = 51.5110*pi/180;
dec2 = 27.9911*pi/180;
dec3 = 14.6609*pi/180;

st1 = 0.0*pi/180;
st2 = 0.250684*pi/180;
st3 = 0.501369*pi/180;

Rphi = Re/sqrt(1-(2*flat - flat^2)*sin(phi)^2);
Rc = Rphi + H;
Rs = (1-flat)^2*Rphi + H;

Rx = Rc*cos(phi);
Ry = Rx;
Rz = Rs*sin(phi);

R = [Rx*cos(st1),Ry*sin(st1),Rz;Rx*cos(st2),Ry*sin(st2),Rz;Rx*cos(st3),...
    Ry*sin(st3),Rz];


L1X = cos(ra1)*cos(dec1);
L1Y = sin(ra1)*cos(dec1);
L1Z = sin(dec1);
L2X = cos(ra2)*cos(dec2);
L2Y = sin(ra2)*cos(dec2);
L2Z = sin(dec2);
L3X = cos(ra3)*cos(dec3);
L3Y = sin(ra3)*cos(dec3); 
L3Z = sin(dec3);

tau1 = t1 - t2;
tau3 = t3 - t2;
tau = t3 - t1;

rohat1 = [cos(ra1)*cos(dec1),sin(ra1)*cos(dec1),sin(dec1)];
rohat2 = [cos(ra2)*cos(dec2),sin(ra2)*cos(dec2),sin(dec2)];
rohat3 = [cos(ra3)*cos(dec3),sin(ra3)*cos(dec3),sin(dec3)];

p = [cross(rohat2,rohat3);cross(rohat1,rohat3);cross(rohat1,rohat2)];

D = zeros(3,3);
din0 = dot(rohat1,cross(rohat2,rohat3));

for x = 1:3
    for y = 1:3
        D(x,y) = dot(R(x,1:3),p(y,1:3));
    end
end

A = 1/din0*(-D(1,2).*tau3/tau+D(2,2)+D(3,2).*tau1/tau);
B = 1/(6*din0)*(D(1,2).*(tau3^2-tau^2)*tau3/tau+D(3,2).*(tau^2-tau1^2)*...
    tau1/tau);file://///iastate.edu/cyfiles/inorris/Desktop/html/astro451hw4set5.html
E = dot(R(2,1:3),rohat2);

a = -(A.^2+2.*A.*E+(R(2,1)^2+R(2,2)^2+R(2,3)^2));
b = -2*mu.*B.*(A+E);
c = -mu^2.*B.^2;

syms q
rsub2 = solve((q.^8+a*q.^6+b*q.^3+c)==0,q);
v = 1;
r2a = double(rsub2);

for l = 1:8
    if r2a(l)>Re
        r2 = r2a(l);
    end
end

ro1 = 1/din0*((6*(D(3,1)*tau1/tau3+D(2,1)*tau/tau3)*r2^3 + mu*D(3,1)...
    *(tau^2-tau1^2)*tau1/tau3)/(6*r2^3+mu*(tau^2-tau3^2))-D(1,1));
ro2 = A + mu*B/r2^3;
ro3 = 1/din0*((6*(D(1,3)*tau3/tau1+D(2,3)*tau/tau1)*r2^3 + mu*D(1,3)...
    *(tau^2-tau3^2)*tau3/tau1)/(6*r2^3+mu*(tau^2-tau1^2))-D(3,3));
ro3 = 963.2744;

rvec1 = [R(1,1),R(1,2),R(1,3)] + ro1*rohat1;
rvec2 = [R(2,1),R(2,2),R(2,3)] + ro2*rohat2;
rvec3 = [R(3,1),R(3,2),R(3,3)] + ro3*rohat3;

f1 = 1 - .5*mu/r2.^3*tau1^2;
f3 = 1 - .5*mu/r2^3*tau3^2;
g1 = tau1 - 1/6*mu/r2^3*tau1^2;
g3 = tau3 - 1/6*mu/r2^3*tau3^2;

v2 = 1/(f1*g3-f3*g1)*(-f3*rvec1+f1*rvec3);
v2mag = norm(v2);

display(r2)
display(v2mag)






