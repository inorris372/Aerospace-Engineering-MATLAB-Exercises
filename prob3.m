%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Astrodynamics 2 # 3         %
%      Ian Norris                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial Conditions

t1 = 0;
t2 = 3;
tspan = linspace(t1,t2,.01);
ta = 20;
x = 0;
l = .1;
g = 9.81;
thetad = 0;
m1 = 4;
m2 =3;

Ic(1) = ta;
Ic(2) = l;
Ic(3) = g;
Ic(4) = thetad;
% syms theta(t) x(t)
% V = odeToVectorField(diff(theta,2) == -(((-g*l*sin(theta)-l^2*diff...
%     (theta,2))/(l*cos(theta)))*cos(theta)+g/l*sin(theta)))
% B = odeToVectorField(diff(x,2) == (m2*l*(-(diff(x,2)*cos(theta)+g/l*...
%     sin(theta)))*cos(theta))/(m1+m2))
%-diff(theta)^2*sin(theta)
options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10);
[theta, x] = ode45(@motion, tspan, [Ic(1) Ic(2) Ic(3) Ic(4)], options);

