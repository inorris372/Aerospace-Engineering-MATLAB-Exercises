%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Lunar Trajectory Project - Ian Norris    %
%          Patched Conics Trajectory          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%It should be noted that for the first part of this project all
%values will be acquired using canonical units.
%1 TU = 806.80415 s
%1 DU = 6378.1km
 
%Constants
clc
clear home
 
RE = 6378.1;  %Radius of Earth (km)
RM = 1738;  %Radius of Moon (km)
Rsoi = 66183;  %Radius of the Sphere of Influence of the Moon (km) 
RsoiA = Rsoi/RE; %(DU)
D = 384400;  %Distance between the Earth and Moon (km)
Da = D/RE; %(DU)
MuE = 3.9860044*10^5; MuE1 = 1;  %Earth Gravitational Parameter (km^3/s^2)
wM = 2.6491065*10^-6;  %Radial Velocity of the Moon (rad/s)
muM = 4.9028*10^3;  %Moon Gravitational Parameter (km^3/s^2)
LOE = 185;  %Parking orbit (km)
Lalt = 120.1; %Periselenium Altitude (km)
phi0 = 0; %Starting value for Flight Path Angle (DEG)
lambda1 = pi/6; %Starting value for Lambda1 (DEG)
r0 = (RE + LOE)/RE;  %Initial radius prior to departure (DU)
r0a = r0*RE; %(km)
v0 = 1.387177165; %Initial value for Transfer Orbit Velocity
v0a = v0*(RE/806.80415); %Initial Transfer Orbit Velocity(km/s)
vM = wM*D; %Velocity of the Moon relative to the Earth (km/s)
 
 
% phi0 = zeros(10000);  
% lambda1 = zeros(10000);  
% v0 = zeros(10000); 
% Maxv0 = sqrt(2*MuE1/r0);  %Upper Limit for v0 (DU/TU)
% MaxLambda1 = pi;  %Upper Limit for Lambda1 (rad)
% Maxphi0 = pi/2;  %Upper Limit for Flight Path Angle (rad)
% v0 = 1.3778;  %Zeroed Initial Velocity prior to departure (DU/TU)
% for x = 1:10000
%     lambda1(x) = (x/10000)*MaxLambda1;
% end
% for y = 1:10000
%     phi0(x) = (y/10000)*Maxphi0;
% end
% for z = 1:10000
%     v0(z) = (z/10000)*Maxv0;
% end
% 
% 
% surf(v0,r0,lambda1)
% %     for 0:1:MaxLambda1
% %         for 0:1:Maxphi0
            
    
%Variables
 
Hm = r0*v0*cos(phi0); %(DU^3/TU^2)
Hma = r0a*v0a*cos(phi0); %(km^2/s^2)
%Angular Momentum as a function of Flight Path Angle
Em0 = (v0^2)/2-MuE1/r0; %(DU^2/TU^2)
Em0a = (v0a^2)/2-MuE/r0a; %(km^2/s^2)
%Initial Specific Mechanical Energy of the Spacecraft in its
%parking orbit
r1 = sqrt(Da^2+RsoiA^2-2*Da*RsoiA*cos(lambda1)); %(DU)
r1a = sqrt(D^2+Rsoi^2-2*D*Rsoi*cos(lambda1)); %(km)
%Radius from Earth Center to Moon Sphere of Influence as
%a function of lambda1
v1 = sqrt(2*(Em0+MuE1/r1)); %(DU/TU)
v1a = sqrt(2*(Em0a+MuE/r1a)); %(km/s)
%Velocity of Spacecraft from Earth parking orbit to Moon Sphere
%of Influence
phi1 = acos(Hm/(r1*v1));
%Flight path angle of the Spacecraft while on its trajectory to
%the Moon's Sphere of Influence
gamma1 = asin(Rsoi/r1a*sin(lambda1));
%Angle between the r1 vector and Distance vector between the 
%Earth and the Moon
p1 = Hm^2/MuE1;
p1a = Hma^2/MuE;
%Semi-Latus Rectum of the Spacecraft transfer orbit
a1 = -MuE1/(2*Em0); %(DU)
a1a = -MuE/(2*Em0a); %(km)
%Semi-Major Axis length of the transfer orbit
e1 = sqrt(1-p1a/a1a);
%Eccentricity of the transfer orbit
V0 = acos((p1-r0)/(r0*e1));
if V0<1*10^-5
    V0=0;
end
%True Anomaly of parking orbit
V1 = acos((p1-r1)/(r1*e1));
%True Anomaly of Transfer orbit to Moon SOI
E0 = acos((e1+cos(V0))/(1+e1*cos(V0)));
%Eccentric Anomaly of Parking Orbit
E1 = acos((e1+cos(V1))/(1+e1*cos(V1)));
%Eccentric Anomaly of Transfer Orbit
TOF = sqrt((a1a^3)/MuE)*((E1-e1*sin(E1))-(E0-e1*sin(E0)));
%Time of Flight to Moon SOI
B = TOF * wM;
%Radians the Moon has traveled
 
 
%Now Using SI Units
 
gamma0 = V1-V0-gamma1-B;
%Phase Angle at Departure
r2 = Rsoi;
%Radius 2 is equivalent to the radius at the Moon's Sphere of
%Influence
v2 = sqrt(v1a^2+vM^2-2*v1a*vM*cos(phi1-gamma1));
%Velocity just inside the Moon's Sphere of Influence
epsilon2 = asin(vM/v2*cos(lambda1)-v1a/v2*cos(lambda1+gamma1-phi1));
%Epsilon2 is the angle that defines the direction of the initial
%selenocentric velocity relative to the Moon's Center.
Em2 = v2^2/2-muM/r2;
%Specific Mechanical Energy inside the Moon's SOI.
Hm3 = r2*v2*sin(epsilon2);
%Angular Momentum of Spacecraft inside the Moon's SOI.
p3 = Hm3^2/muM;
%Semi-Latus Rectum of the Spacecraft's orbit about the moon
e3 = sqrt(1+2*Em2*Hm3^2/muM^2);
%Eccentricity of the Spacecraft's orbit about the moon
a3 = p3/(1-e3^2);
%Semi-major Axis of Spacecraft's orbit about the moon
r3 = p3/(1+e3);
%Final distance of the Spacecraft from the moon's center
v3 = sqrt(2*(Em2+muM/r3));
%Final Velocity of the Spacecraft at Pereselenium
phi2 = acos(Hm3/(r3*v3));
%Flight path angle of the Spacecraft in Moon orbit
if phi2<1*10^-5
    phi2 = 0;
end
V3 = acos((p3-r2)/(r2*e3));
%True Anomaly of Spacecraft just inside Moon SOI
V4 = acos((p3-r3)/(r3*e3));
%True Anomaly at Periselenium of Moon
F1 = 2*atanh(sqrt((e3-1)/(e3+1))*tan(V3/2));
%Eccentric Anomaly of Spacecraft just inside Moon SOI
F2 = 2*atanh(sqrt((e3-1)/(e3+1))*tan(V4/2));
%Eccentric Anomaly of Spacecraft at Periselenium of Moon
TOF2 = sqrt(((-a3)^3)/muM)*((e3*sinh(F1)-F1)-(e3*sinh(F2)-F2));
%Time of Flight from Moon Sphere of Influence to Periselenium
 
fprintf('Astrodynamics Lunar Trajectory Results\n')
fprintf('v0 = %g (DU/TU), phi0 = %g (DU^2/TU)\n',v0,phi0*180/pi)
fprintf('lambda1 = %g (DEG), Em0 = %g (DU^2/TU^2)\n',lambda1*180/pi,Em0)
fprintf('V0 = %g (DEG), V1 = %g (DEG)\n',V0*180/pi,V1*180/pi)
fprintf('B = %g (DEG), gamma0 = %g (DEG)\n',B*180/pi,gamma0*180/pi)
fprintf('r1 = %g (DU), v1 = %g (DU/TU)\n',r1,v1)
fprintf('phi1 = %g (DEG), gamma1 = %g (DEG)\n',phi1*180/pi,gamma1*180/pi)
fprintf('a1 = %g (DU), e1 = %g\n',a1,e1)
fprintf('v2 = %g (km/s), phi2 = %g (DEG)\n',v2,phi2*180/pi)
fprintf('epsilon2 = %g (DEG), Em2 = %g (km^2/s^2)\n',epsilon2*180/pi,Em2)
fprintf('a3 = %g (km), e2 = %g\n',a3,e3)
fprintf('r3 = %g (km), v3 = %g (km/s)\n',r3,v3)
fprintf('Time of Flight from Parking Orbit to Moon SOI\n')
fprintf('TOF = %g (Hours)\n',TOF/3600)
fprintf('Time of Flight from Moon SOI to Periselenium\n')
fprintf('TOF2 = %g (Hours)\n',TOF2/3600)
