%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Ian Norris           %
%           Midterm Exam          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
%% Problem 1
% General Variable Initialization
Re = 6378.1;             %[km] Radius of Earth
mu = 3.986e5;            %[km^3/s^2] Gravitational Parameter
vnom = pi/2;             %[rad] True Anomaly Angle
vnom1 = 0;               %[rad] True Anomaly Angle point 1
vnom2 = pi;              %[rad] True Anomaly Angle point 2
tol = 10e-20;            % Error Tolerance
told = 1;                % Counter

% Hyperbolic Orbit Variable Initialization
ahyp = -3e4;             %[km]  Semi-Major Axis
ehyp = sqrt(2);          % Unitless Eccentricity
del = asin(1/ehyp)*2;    %[rad] Turning Angle
r0phyp = ahyp*(1-ehyp);  %[km] Radius at Perigee
phyp = ahyp*(1-ehyp^2);  %[km] Semi-latus Rectum


% Elliptic Orbit Variable Initialization
ael = 15000;             %[km]  Semi-Major Axis
eel = .5;                % Unitless Eccentricity
rael = eel*ael+ael;      %[km]  Radius of Apogee
E1 = acos((eel + cos(vnom1))/(1+eel*cos(vnom1)));  % Eccentric Anomaly
E2 = acos((eel + cos(vnom2))/(1+eel*cos(vnom2)));  % Eccentric Anomaly
tdesire = sqrt(ael^3/mu)*((E2-eel*sin(E2))-(E1-eel*sin(E1)));
%Transfer  Orbit time from point 1 to point 2 that is desired

r1 = [0 phyp 0];         %[km] Radius 1
r2 = [-rael 0 0];        %[km] Radius 2
%Gauss Variables
A = sqrt(norm(r1)*norm(r2))*sin(vnom)/sqrt(1-cos(vnom));
z = linspace(-2*pi,4*pi^2,1000000);
C = .5-z/factorial(4)+z.^2/factorial(6)-z.^3/factorial(8)+z.^4/...
    factorial(10);
S = 1/3-z/factorial(5)+z.^2/factorial(7)-z.^3/factorial(9)+z.^4/...
    factorial(11);
Y = norm(r1) + norm(r2) - A*(1-z.*S)./sqrt(C);
x = sqrt(Y./C);
t = (x.^3.*S+A.*sqrt(Y))./sqrt(mu);    %Actual Time of Flight (1-2)
error = abs(tdesire-t);       %Difference between Actual and Desired Time
[b,ersmall] = min(error(:));     %Minimum Difference in Time
zwant = z(ersmall);           %Z that is found using Minimum Difference
%zwant = .6568;
%New Gauss Variables
Cv = .5-zwant/factorial(4)+zwant^2/factorial(6)-zwant^3/factorial(8)+...
    zwant^4/factorial(10);
Sv = 1/3-zwant/factorial(5)+zwant^2/factorial(7)-zwant^3/factorial(9)+...
    zwant^4/factorial(11);
Yv = norm(r1) + norm(r2) - A*(1-zwant*Sv)/sqrt(Cv);
%Lagrange F & G Coefficients
f = 1 - Yv/norm(r1);
g = A*sqrt(Yv/mu);
gdot = 1 - Yv/norm(r2);
v1 = (r2-f.*r1)./g;
v2 = (gdot.*r2-r1)./g;
% Hyperbolic Velocity
vH = sqrt(2*mu./r1-mu/ahyp);
% Elliptic Velocity
vE = sqrt(2*mu./r2-mu/rael);
e = 1./mu.*((v1.^2-mu./norm(r1)).*r1 - dot(r1,v1).*v1); 
%Eccentricity Vector
a = norm(r2)/(norm(e)+1);     %Semi-Major Axis
delv1 = v1-vH;                    %Delta V at point 1
delv2 = vE-v2;                    %Delta V at point 2
%Output all Variable Values
display(v1)
display(v2)
display(e)
display(a)
display(delv1)
display(delv2)


%% Problem 2

% Circular Parking Orbit
Re = 6378.1;             %[km] Radius of Earth
mu = 3.986e5;            %[km^3/s^2] Gravitational Parameter
r1 = 6700;               %[km] Parking Orbit Radius
a1 = r1;                 %[km] Semi-Major Axis
ee1 = 0;                 % Eccentricity
i1 = 28.5*pi/180;        %[rad] Initial Inclination Angle
vc1 = sqrt(mu*(2/r1-1/a1));  %[km/s] Initial Velocity

% Geostationary Orbit
r2 = 42164;              %[km] Radius of Geostationary Orbit
i2 = 0;                  %[rad] GTO Inclination Angle
a2 = r2;                 %[km] GTO Semi-Major Axis
ee2 = 0;                 % Eccentricity
vc2 = sqrt(mu*(2/r2-1/a2));  %[km/s] Final Velocity
% General Variables

theta = abs(i1-i2);      %[rad] Change in Inclination Angle
%Initialize Storage Variables
cstore = 100000000;      % Storage variable for 1st Cost Function
cstore2 = cstore;        % Storage variable for 2nd Cost Function

% Optimal Hohmann Transfer With Plane Change
rptr = r1;               %[km] Radius at Perigee of Transfer Orbit
ratr = r2;               %[km] Radius at Apogee of Transfer Orbit
at = (rptr + ratr)/2;    %[km] Semi-Major Axis of Transfer Orbit
et = (ratr - at)/at;     % Eccentricity
vptr = sqrt(2*mu/rptr-mu/at);  %[km/s] Velocity at Perigee of Transfer
vatr = sqrt(2*mu/ratr-mu/at);  %[km/s] Velocity at Apogee of Transfer
delv1H = vptr - vc1;     %[km/s] Delta V needed at Perigee of Transfer
delv2H = vatr - vc2;     %[km/s] Delta V needed at Apogee of Transfer

for alpha1 = i2:.0001:i1            %Use Loop to vary Alpha 1&2
    alpha2 = theta - alpha1;        %Initialize Alpha 1&2 values
    %Transfer Orbit Equations
    delV1 = sqrt(vptr^2 + vc1^2 - 2*vptr*vc1*cos(alpha1));       %[km/s] 
    delV2 = sqrt(vc2^2 + vatr^2 - 2*vc2*vatr*cos(theta-alpha1)); %[km/s]
    delVT = abs(delV1) + abs(delV2); %[km/s] Total Delta V needed
    v1 = pi;                         %[rad] True Anomaly Angle
    ToF1 = pi*sqrt(at^3/mu)/3600;    %[Hours] Time of Flight
    C1 = delVT + 2*ToF1;             %[Unitless] Cost Function
    if C1 < cstore
        alstoreA = alpha1;
        alstoreB = alpha2;
        cstore = C1;
        nustore = v1;
    end
    
end

% One Tangent Burn
for alpha1b = i2:.01:i1             %Use Loop to vary Alpha 1&2
    alpha2b = theta - alpha1b;      %Initialize Alpha 1&2 values
    for at2 = at:1:50000            %Vary Semi-Major Axis via second loop
        %Transfer Orbit Equations
        et2 = 1 - ratr/at2;         %Eccentricity
        nu2 = acos((at2*(1-et2^2)-ratr)/(et2*ratr));
        %[rad] Initial True Anomaly Angle
        phi = atan(et2*sin(nu2)/(1+et2*cos(nu2)));
        %[rad] Initial Flight Path Angle
        vatr2 = sqrt(2*mu/ratr-mu/at2);  %[km/s] Velocity at Apogee
        delVAT = vatr2 - vc1;            %[km/s] First Delta V
        %delVAT = sqrt(vptr2^2 + vc1^2 - 2*vptr2*vc1*cos(alpha1b));
        vn2 = acos((at2*(1-et2^2)-rptr)/(et2*rptr));
        %[rad] Final True Anomaly Angle
        vptr2 = sqrt(2*mu/rptr-mu/at2);  %[km/s] Velocity at Perigee
        delVBT = sqrt(vc2^2 + vptr2^2 - 2*vc2*vptr2*cos(theta-alpha1b));
        %[km/s] Final Delta V
        delVT2 = abs(delVAT) + abs(delVBT);   %[km/s] Total Delta V
        %[rad] Eccentric Anomaly
        E2 = atan(sqrt(1-et2^2)*sin(vn2)/(et2+cos(vn2)));
        ToF2 = (E2-et2*sin(E2))*sqrt(at2^3/mu)/3600;  
        %[hours] Final Time of Flight 
    
        %Cost Function
        C2 = delVT2 + 2*ToF2;
        if C2 < cstore2
            % Storage Variables
            alstore1 = alpha1b;
            alstore2 = alpha2b;
            cstore2 = C2;
            nu2store = nu2;
            atstore = at2;
        end
    end
end
%Plot and Display User Required Variables
display(ToF1)
display(alstore1)
display(alstore1)
display(nustore)
display(cstore)
display(ToF2)
display(alstore1)
display(alstore2)
display(nu2store)
display(cstore2)
figure(1)

figure(2)

%% Problem 3

%First initialize all necessary variables
r0v = [20000 -105000 -19000]; % Initial radius (kilometers)
v0v = [.9 -3.4 -1.5]; % Initial velocity (kilometers/sec.)
r0 = norm(r0v); %km   (Initial radius)
v0 = norm(v0v); %(km/s)    (Initial velocity)

mu = 3.986e5; % (Earth's Gravitational Parameter) units = (km^3*s^-2)
Em = v0^2/2-mu/r0; % Specific Mechanical Energy
a = 1/(2/r0-v0^2/mu);  % Semi-Major Axis
Hm = cross(r0v,v0v);   % Angular Momentum Vector
ee = sqrt(1+2*Em*norm(Hm)^2/mu^2);  % Eccentricity
nu0 = acos((a*(1-ee^2)-r0)/(ee*r0));  %Initial True Anomaly Angle
T0 = 0;  %Initial Time (seconds)
t = 2*60*60;  % Final Time, so 2 * 60 minutes * 60 seconds
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
c(n) = .5 - z(n)/factorial(4)+z(n)^2/factorial(6)-z(n)^3/factorial(8)...
  + z(n)^4/factorial(10)-z(n)^5/factorial(12);
s(n) = 1/6 - z(n)/factorial(5)+z(n)^2/factorial(7)-z(n)^3/factorial(9)...
  + z(n)^4/factorial(11)-z(n)^5/factorial(13);

f(n)=((1-r0/a)*s(n)*x(n)^3+dot(r0v,v0v)/sqrt(mu)*c(n)*x(n)^2+r0*x(n)-sqrt(mu)*t);
%Initial Function of Universal Variable
fp(n)=c(n)*x(n)^2+dot(r0v,v0v)/sqrt(mu)*(1-s(n)*z(n))*x(n)+r0*(1-c(n)*z(n));
%First derivative of Universal Variable Function
fpp(n)=(1-r0/a)*(1-s(n)*z(n))*x(n)+dot(r0v,v0v)/sqrt(mu)*(1-c(n)*z(n));
%Second derivative of Universal Variable Function
deln = 2*sqrt(abs(4*fp(n)^2-5*f(n)*fpp(n)));
delX(n)= 5*f(n)/(fp(n)+ sign(fp(n))*deln);
%fp(n)/abs(fp(n))*deln);
%Change in Universal Variable
x(n+1) = x(n)-delX(n);
%New value of Universal Variable
%Lower If statement shuts off while loop upon reaching epsilon
if(abs(delX(n)^2/a)<epsilon)
    flag = false;
end
%Counter is raised below each iteration
n = n+1;
end

%Plug in Universal Variable into f & g equations to find new r and v values

fla = 1 - x(n-1)^2/r0*c(n-1);
gla = t - x(n-1)^3/sqrt(mu)*s(n-1);
r = fla*r0v+gla*v0v;
fdot = sqrt(mu)/(norm(r)*r0)*(s(n-1)*z(n-1)-1)*x(n-1);
gdot = 1 - x(n-1)^2/norm(r)*c(n-1);
v = fdot*r0v+gdot*v0v;

Emcheck = norm(v)^2/2-mu/norm(r); % Specific Mechanical Energy

display(r)
display(v)
display(Emcheck)








