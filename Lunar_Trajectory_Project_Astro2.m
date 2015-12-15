%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Lunar Trajectory Project - Ian Norris    %
%          Patched Conics Trajectory          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 

%Constants
clc
clear home
global mu wM
RE = 6378.1;  %Radius of Earth (km)
RM = 1738;  %Radius of Moon (km)
D = 384400;  %Distance between the Earth and Moon (km)
MuE = 3.9860044*10^5;  %Earth Gravitational Parameter (km^3/s^2)
wM = 2.6491065*10^-6;  %Radial Velocity of the Moon (rad/s)
muM = 4.9028*10^3;  %Moon Gravitational Parameter (km^3/s^2)
LEO = 222;  %Parking orbit (km)
vp = sqrt(MuE*(1/(RE+LEO)));
phi0 = 0; %Starting value for Flight Path Angle (DEG)
vM = wM*D; %Velocity of the Moon relative to the Earth (km/s)
mu = muM/MuE;
x1 = mu;
x2 = -1 + mu;
y1 = 0;
y2 = 0;
tof = 10;
dt = .003;
t2 = -dt;
vecT = linspace(0,tof,dt);
vecT2(1:size(vecT)) = 1e-10; 
%tof = 9*24*60*60;
%[ R ] = rstTOijk(x)
gamma = 90;
%75:150
mag = 10.98;
%Delta V = 10.6:10.98
x(1) = x1*D-(RE+LEO)*cosd(gamma);
dxbydt(1) = -(vp+mag)*sind(gamma);
y(1) = (RE+LEO)*sind(gamma);
dybydt(1) = (vp+mag)*cosd(gamma);
options = odeset('RelTol',1e-10,'AbsTol',[1e-10 1e-10 1e-10 1e-10]);
npt = 0;
tf = pi;
t1 = t2;
t2 = t1 + dt;
[trk,ysol] = ode45(@EquationsOfMotion,[0.0:dt:1000],[x(1) y(1) dxbydt(1)...
    dybydt(1)],options);
npt = npt + 1;
xplot(npt) = ysol(length(trk), 1);
yplot(npt) = ysol(length(trk), 2);

%% Third Lagrange Point
 
%0 = x + (1-mu)/(mu+x)^2 + mu/(x-1+mu)^2;
% xup = 2;
% xdown = -2;
% % L3 libration point
% 
% ilp = 3;
% xr1 = -2;
% xr2 = +2;
% rtol = 1.0e-8;
% [xl3, froot] = brent ('clpfunc', xr1, xr2, rtol);
% yl3 = 0;
% r1sqr = (xl3 - xm1)^2 + yl3^2;
% r2sqr = (xl3 - xm2)^2 + yl3^2;
% e3 = -0.5 * (xl3^2 + yl3^2) - (1 - mu) / sqrt(r1sqr) - mu / sqrt(r2sqr);

%% Plot Circular Restricted Three Body Problem

figure(1)

[r1,s1] = getCircle(x2*D,y1*D,RM);
[r2,s2] = getCircle(x2*D,y1*D,RM+100);
[r3,s3] = getCircle(x2*D,y1*D,RM+200);
[r4,s4] = getCircle(x1*D,y2*D,RE);
[r5,s5] = getCircle(x1*D,y2*D,RE+150);
[r6,s6] = getCircle(x1*D,y2*D,RE+200);
[r7,s7] = getCircle(x1*D,y2*D,RE+250);
[r8,s8] = getCircle(x1*D,y2*D,D);

H1 = patch(r1,s1,1); %Make Moon patch object
set(H1,'FaceColor',[1,1,1]*0.75) % Set the Moon Color to Grey
H2 = patch(r4,s4,1); %Make Earth patch object 
set(H2,'FaceColor',[0,0,1]*0.75) % Set the Earth Color to Blue

plot(r1,s1,'-b','linewidth',1)
hold on
plot(r2,s2,'-y','linewidth',2)
hold on
plot(r3,s3,'-g','linewidth',2)
hold on
plot(r4,s4,'-b','linewidth',1)
hold on
plot(r5,s5,'-r','linewidth',2)
hold on
plot(r6,s6,'-y','linewidth',2)
hold on
plot(r7,s7,'-g','linewidth',2)
hold on
plot(r8,s8,'-r','linewidth',1)
hold on
plot(xplot,yplot,'-b')
grid on
axis square
% axis ([-.7e5 10e5 -.1e5 1.5e5]) 









