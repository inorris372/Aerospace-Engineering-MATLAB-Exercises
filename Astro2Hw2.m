%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Ian Norris                %
%    AerE 451 Hw Set 2         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1

% Find the local sidereal time (in degrees) at east longitude 139.8(deg) on 
% March 3, 2004 at 4:30:00 UT.[Hint: answer lies between 0 and 20(deg)]

%First calculate the Julian Date.

year = 2004;
month = 3;
day = 3;
lambdaE = 139.8; % Degrees
UT = 4.5;

Jnot=367*year-floor(7/4*(year+floor((month+9)/12)))+floor(275*(month/9))...
    + day + 1721013.5;

Tnot = (Jnot - 2451545.0)/36525;

ThetaGnot = 100.4606184 + 36000.77004*Tnot + 0.000387933*Tnot^2 ...
    - 2.583e-8*Tnot^3;
ThetaG = ThetaGnot + 360.98564724*UT/24;
Theta = ThetaG + lambdaE;
ThetaP = floor(Theta/360);
ThetaDone = (Theta/360 - ThetaP)*360;
display(ThetaDone)

%% Problem 2

% An Earth satellite has an initial radius of r0 = 10000 km and an initial
% true anomaly angle of nu0 = 30(deg) at a speed of v0 = 10 km/s. Use the 
% Universal Variable approach to Kepler’s problem to solve for the true 
% anomaly angle nu1 one hour after the initial observations. Record all 
% intermediate results and show all work.

%First initialize all necessary variables
r0 = 1000; %km   (Initial radius)
nu0 = 30; %(degrees)   (Initial true anomaly angle)
v0 = 10; %(km/s)    (Initial velocity)
mu = 3.986e5; % (Earth's Gravitational Parameter) units = (km^3*s^-2)
Em = v0^2/2-mu/r0; % Specific Mechanical Energy
a = -mu/(2*Em);  % Semi-Major Axis
T0 = 0;  %Initial Time (seconds)
t = 60*60;  % Final Time, so 60 minutes * 60 seconds
deltaT= t - T0; %Change in time that we are measuring (seconds)
n = 1;  % Counter
% Must initialize all array sets
epsilon = 10^-8;   %Error Limit to test for convergence
flag = true;
%z = zeros;
c = zeros;
s = zeros;
f = zeros;
fp = zeros;
fpp = zeros;
delX = zeros;
x = zeros;
z = zeros;

while(flag)
%Initial value of Universal Variable
x(1) = sqrt(mu)*deltaT/abs(a);

%Laguerre Approximation is implemented
z(n) = x(n)^2/a;
c(n) = .5 - z(n)/factorial(4)+z(n)^2/factorial(6)-z(n)^3/factorial(8);
s(n) = 1/6 - z(n)/factorial(5)+z(n)^2/factorial(7)-z(n)^3/factorial(9);

f(n)=(1-r0/a)*s(n)*x(n)^3+r0.*v0./sqrt(mu).*c(n).*x(n)^2+r0*x(n)-sqrt(mu)*t;
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
%Counter is raised below each iteration
n = n+1;
end

%Plug in Universal Variable into f & g equations to find new r and v values
fla = 1 - x(n-1)^2/r0*c(n-1);
gla = t - x(n-1)^3/sqrt(mu)*s(n-1);
r = fla*r0+gla*v0;
fdot = sqrt(mu)/(r*r0)*(s(n-1)*z(n-1)-1)*x(n-1);
gdot = 1 - x(n-1)^2/r*c(n-1);
v = fdot*r0+gdot*v0;

display(r)
display(v)

%% Problem 3

% An Earth satellite moves in the equatorial plane of the ECI frame with an
% initial position and velocity of r0 = 7000.0 Iˆ ? 12 124 Jˆ km, and 
% v0 = 2.6679 Iˆ + 4.6210 Jˆ km/s. Using the Laguerre approximation in the 
% Universal Variable formulation, together with the Lagrange f and g 
% functions, compute the position and velocity vectors of the satellite 60 
% minutes after the initial observations. Show all intermediate work.

%First initialize all necessary variables
r0v = [7000 -12124]; % Initial radius (kilometers)
v0v = [2.6679 4.6210]; % Initial velocity (kilometers/sec.)
r0 = norm(r0v); %km   (Initial radius)
v0 = norm(v0v); %(km/s)    (Initial velocity)
%r0 = 7000;
%v0 = 2.6679;
%r0 = -12124;
%v0 = 4.6210;
mu = 3.986e5; % (Earth's Gravitational Parameter) units = (km^3*s^-2)
Em = v0^2/2-mu/r0; % Specific Mechanical Energy
a = -mu/(2*Em);  % Semi-Major Axis
T0 = 0;  %Initial Time (seconds)
t = 60*60;  % Final Time, so 60 minutes * 60 seconds
deltaT= t - T0; %Change in time that we are measuring (seconds)
n = 1;  % Counter
% Must initialize all array sets
epsilon = 10^-8;   %Error Limit to test for convergence
flag = true;
%z = zeros;
c = zeros;
s = zeros;
f = zeros;
fp = zeros;
fpp = zeros;
delX = zeros;
x = zeros;
z = zeros;

while(flag)
%Initial value of Universal Variable
x(1) = sqrt(mu)*deltaT/abs(a);

%Laguerre Approximation is implemented
z(n) = x(n)^2/a;
c(n) = .5 - z(n)/factorial(4)+z(n)^2/factorial(6)-z(n)^3/factorial(8);
s(n) = 1/6 - z(n)/factorial(5)+z(n)^2/factorial(7)-z(n)^3/factorial(9);

f(n)=(1-r0/a)*s(n)*x(n)^3+dot(r0v,v0v)./sqrt(mu).*c(n).*x(n)^2+r0*x(n)-sqrt(mu)*t;
%Initial Function of Universal Variable
fp(n)=c(n)*x(n)^2+dot(r0v,v0v)./sqrt(mu).*(1-s(n).*z(n)).*x(n)+r0*(1-c(n)*z(n));
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
%Counter is raised below each iteration
n = n+1;
end

%Plug in Universal Variable into f & g equations to find new r and v values
fla = 1 - x(n-1)^2/r0*c(n-1);
gla = t - x(n-1)^3/sqrt(mu)*s(n-1);
r = fla*r0+gla*v0;
fdot = sqrt(mu)/(r*r0)*(s(n-1)*z(n-1)-1)*x(n-1);
gdot = 1 - x(n-1)^2/r*c(n-1);
v = fdot*r0+gdot*v0;

display(r)
display(v)





