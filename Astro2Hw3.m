%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ian Norris                  %
%      Astrodynamics 451           %
%      Homework Set 3              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1

% 1. An Earth satellite has a perigee radius of 7000 km and an apogee 
% radius of 10 000 km. What is the value of the true anomaly angle V swept 
% out during propagation from 30 minutes after perigee passage to 90 
% minutes after perigee? Use Universal Variables and Laguerre solution.
Rp = 7000; %km   (Perigee radius)
Ra = 10000; %km   (Apogee radius)
a = (Rp + Ra)/2;   % Semi-Major Axis   (km)
e = (Ra/a)-1;  % Eccentricity   (unitless)
mu = 3.986e5; % (Earth's Gravitational Parameter) units = (km^3*s^-2)
Em = mu/(2*a);% Specific Mechanical Energy
Vp = sqrt(2*(Em + mu/Rp)); % Velocity at Perigee  (km/s)
Va = Vp*Rp/Ra;  % Velocity at Apogee (km/s)
T0 = [30*60 90*60];  %Time that has passed (seconds)

r0 = Rp; %km   (Initial radius)
nu0 = 0; %(degrees)   (Initial true anomaly angle)
v0 = Vp; %(km/s)    (Initial velocity)

% Must initialize all array sets
epsilon = 10^-8;   %Error Limit to test for convergence
flag = true;
c = zeros;
s = zeros;
f = zeros;
fp = zeros;
fpp = zeros;
delX = zeros;
x = zeros;
z = zeros;
deltaT = [0 0];
nuF = [0 0];

for b = 1:2

deltaT(b)= T0(b); %Change in time that we are measuring (seconds)
n = 0;  % Counter
while(flag)
%Counter is raised below each iteration
n = n+1;
%Initial value of Universal Variable
x(1) = sqrt(mu)*deltaT(b)/abs(a);

%Laguerre Approximation is implemented
z(n) = x(n)^2/a;
c(n) = .5 - z(n)/factorial(4)+z(n)^2/factorial(6)-z(n)^3/factorial(8)+...
    z(n)^4/factorial(10)-z(n)^5/factorial(12);
s(n) = 1/6 - z(n)/factorial(5)+z(n)^2/factorial(7)-z(n)^3/factorial(9)+...
    z(n)^4/factorial(11)-z(n)^5/factorial(13);

f(n)=(1-r0/a)*s(n)*x(n)^3+r0.*v0./sqrt(mu).*c(n).*x(n)^2+r0*x(n)-sqrt(mu)*T0(b);
%Initial Function of Universal Variable
fp(n)=c(n)*x(n)^2+r0.*v0./sqrt(mu).*(1-s(n).*z(n)).*x(n)+r0*(1-c(n)*z(n));
%First derivative of Universal Variable Function
fpp(n)=(1-r0/a)*(1-s(n)*z(n))*x(n)+r0.*v0./sqrt(mu).*(1-c(n).*z(n));
%Second derivative of Universal Variable Function
deln = 2*sqrt(4*(fp(n))^2-5*f(n)*fpp(n));
delX(n)=5*f(n)/(fp(n)+sign(fp(n))*deln);
%Change in Universal Variable
x(n+1) = x(n)-delX(n);
%New value of Universal Variable
%Lower If statement shuts off while loop upon reaching epsilon
if((delX(n)^2/a)<epsilon)
    flag = false;
end
%Plug in Universal Variable into f & g equations to find new r value

fla = 1 - x(n)^(2)/r0*c(n);
gla = T0(b) - x(n)^(3)/sqrt(mu)*s(n);
end

r = fla*r0+gla*v0;

%Calculate Anomaly Angle for each time duration
nuF(b) = acos((a*(1-e^2)-r)/(r*e));
end
%Calculate total Anomaly Angle from difference between two Anomaly Angles
%that were previously calculated.
nuDone = (nuF(2)-nuF(1))*180/pi;
display(nuDone)

%% Problem 2

% 2. A tracking station at geodetic latitude -20(deg) and elevation 500 m 
% makes the following observations of an Earth satellite:

% L.Time (min) S. Time (_) Azimuth (_) Elevation (_) Range (km)
% 0            60.0        165.932     8.81952       1212.48
% 2            60.5014     145.970     44.2734       410.596
% 4            61.0027     2.40973     20.7594       726.464
% (1)	 Use Reference Ellipsoid. Use Gibbs' method to find the vectors r 
% and v at the 2-minute point, and from these values calculate the orbital 
% elements. (Hint: |v| = 7:7 km/s.)
Req = 6378.145;
phi = -20*pi/180;
eee = .08182;
Elevation = [8.81952 44.2734 20.7594];
Azimuth = [165.932 145.970 2.40973];
SidT = [60.0 60.5014 61.0027];
Rmag = [1212.48 410.596 726.464];
del = zeros(3);
h = zeros(3);
alpha = zeros(3);
x = (Req/sqrt(1-eee^2*sin(phi)^2)+.5)*cos(phi);
z = (Req*(1-eee^2)/sqrt(1-eee^2*sin(phi)^2)+.5)*sin(phi);

for u = 1:3
    del(u) = asin(cos(phi)*cos(Azimuth(u))*cos(Elevation(u))+sin(phi)*...
        sin(Elevation(u)));
    h(u) = 2*pi-acos((cos(phi)*sin(Elevation(u))-sin(phi)*cos(Azimuth(u))*...
    cos(Elevation(u)))/cos(del(u)));
    alpha(u) = SidT(u) - h(u);
    Lx = cos(alpha(u))*cos(del(u));
    Ly = sin(alpha(u))*cos(del(u));
    Lz = sind(del(u));
    rohat = [Lx Ly Lz];
    Rf = [x*cos(SidT(u)*pi/180) x*sin(SidT(u)*pi/180) z];
    R = Rf+rohat.*Rmag(u);
    if u==1
        R1 = R;
    elseif u==2
        R2 = R;
    else
        R3 = R;
    end
end

mu = 398600;               %[km^3/s^2] Earth’s gravitational parameter
r1 = norm(R1); r2 = norm(R2); r3 = norm(R3);
% Vector cross products
R12 = cross(R1,R2);
R23 = cross(R2,R3);
R31 = cross(R3,R1);

Nv  = r1*R23 + r2*R31 + r3*R12;
Dv  = R23 + R31 + R12;
Sv  = (r2-r3)*R1 + (r3-r1)*R2 + (r1-r2)*R3;

N   = norm(Nv);
D   = norm(Dv);
% Velocity vector
V2 = (mu/(N*D))^0.5*(cross(Dv,R2)/r2 + Sv);
display(V2)
v2 = norm(V2);
Hmom2 = cross(R2,V2);
hmo2 = norm(Hmom2);
p = hmo2^2/mu;
a = mu/(-v2^2+(2*mu/r2));
e = sqrt(1-(p/a));
N=cross([0;0;1],Hmom2);
E=(1/mu)*((v2^2-(mu/r2))*R2-(dot(R2,V2))*V2);
I=[1;0;0];
J=[0;1;0];
K=[0;0;1];
i=acos(dot(Hmom2,K)/hmo2);
ideg=i*180/pi;
Omega=acos(dot(N,I)/norm(N));
Odeg=Omega*180/pi;
w=acos(dot(N,E)/(norm(N)*norm(E)));
wdeg=w*180/pi;
nu0=acos(dot(E,R2)/(norm(E)*r2));
nu0deg=nu0*180/pi;

OE = [a e ideg Odeg wdeg nu0deg];

fprintf('a [km]    e [unitless]   i [deg]  Omega[deg]   omega[deg]   Nu0 [deg] \n');
fprintf('%4.2f     %4.2f         %4.2f     %4.4f    %4.2f    %4.2f  \n',OE);
% Check, Only true anomaly
