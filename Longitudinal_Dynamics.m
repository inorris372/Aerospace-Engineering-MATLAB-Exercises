% input data for finding aircraft longitudinal dynamics modes

g = 32.2;
i = sqrt(-1);

% aircraft data

W = 100000; % lbs
S = 1667; % ft^2
AR = 7;
cbar = 15.4; %ft

ib = 1900;

%operational input

u0 = 733; %ft/s
rho = .000889; % slugs/ft^3
theta0 = 0;

% computed values

CL0 = 2*W/(rho*u0^2*S);
Cw0 = 2*W/(rho*u0^2*S);
Iyy = ib*rho*S*cbar^3/8;
m = W/g;
mu = 2*m/(rho*S*cbar);
tstar = cbar/(2*u0);

% stability derivatives

Cxu = -.0376; Czu = 0; Cmu = 0;
Cxalpha = .14; Czalpha = - 4.9; Cmalpha = -.488;
Cxq = 0; Czq = -4.5; Cmq = -22.9;
Cxadot = 0; Czadot = 4; Cmadot = -4.2;

% compute dimensional values

Xu = rho*u0*S*Cw0*sin(theta0)+1/2*rho*u0*S*Cxu;
Zu = -rho*u0*S*Cw0*cos(theta0)+1/2*rho*u0*S*Czu;
Mu = 1/2*rho*u0*cbar*S*Cmu;

Xw = 1/2*rho*u0*S*Cxalpha; Zw = 1/2*rho*u0*S*Czalpha; Mw = 1/2 *rho*u0*cbar*S*Cmalpha;

Xq = 1/4*rho*u0*cbar*S*Cxq; Zq = 1/4*rho*u0*cbar*S*Czq; Mq = 1/4*rho*u0*cbar^2*S*Cmq;

Xwdot = 1/4*rho*cbar*S*Cxadot; Zwdot = 1/4*rho*cbar*S*Czadot; Mwdot = 1/4*rho*cbar^2*S*Cmadot;

fprintf(' \t\t X \t\t\t\t Z \t\t\t\t M \r\n'); fprintf('u\t %f \t %f \t %f\r\n',Xu,Zu,Mu); fprintf('w\t %f \t\t %f \t %f\r\n',Xw,Zw,Mw);
fprintf('q\t %f \t\t %f \t %f\r\n',Xq,Zq,Mq); fprintf('w.\t %f \t\t %f \t\t %f\r\n',Xwdot,Zwdot,Mwdot);
pause
% put values into A matrix

A = [Xu/m Xw/m 0 -g*cos(theta0); Zu/(m-Zwdot) Zw/(m-Zwdot) (Zq+m*u0)/(m-Zwdot) -W*sin(theta0)/(m-Zwdot);
   1/Iyy*(Mu+Mwdot*Zu/(m-Zwdot)) 1/Iyy*(Mw+Mwdot*Zw/(m-Zwdot)) 1/Iyy*(Mq+Mwdot*(Zq+m*u0)/(m-Zwdot)) -Mwdot*W*sin(theta0)/(Iyy*(m-Zwdot));
   0 0 1 0];

[evec eval] = eig(A)

pause
de = poly(A); sys = tf(1,de); pzmap(sys);
pause
% make eigenvectors dimensionless

evecbar1 = [evec(1,:)/u0; evec(2,:)/u0; evec(3,:)*tstar; evec(4,:)];

evecbar = [evecbar1(:,1)/norm(evecbar1(:,1),inf) evecbar1(:,2)/norm(evecbar1(:,2),inf) evecbar1(:,3)/norm(evecbar1(:,3),inf) evecbar1(:,4)/norm(evecbar1(:,4),inf)];

evecamp = abs(evecbar)
evecphase = angle(evecbar)
pause;
figure;
%calculate flight paths

%short period mode
time = linspace(0,6*pi/(imag(eval(1,1))),51); tmax = length(time);
r1 = evecamp(:,1); phi1 = evecphase(:,1); phi2 = evecphase(:,2);
for k = 1:tmax
   t = time(k);
   xsp(:,k) = r1.*(exp(i*phi1).*exp(eval(1,1).*t)+exp(i*phi2).*exp(eval(2,2).*t));
end


subplot(4,1,1),plot(time,xsp(1,:));ylabel('\^{u}');title('Short-Period Mode');
subplot(4,1,2),plot(time,xsp(2,:));ylabel('\alpha');
subplot(4,1,3),plot(time,xsp(3,:));ylabel('\^{q}');
subplot(4,1,4),plot(time,xsp(4,:));ylabel('\theta');xlabel('time');

figure;
%phugoid mode

time = linspace(0,6*pi/(imag(eval(3,3))),51); tmax = length(time);
r3 = evecamp(:,3); phi1 = evecphase(:,3); phi2 = evecphase(:,4);
for k = 1:tmax
   t = time(k);
   xsp(:,k) = r3.*(exp(i*phi1).*exp(eval(3,3).*t)+exp(i*phi2).*exp(eval(4,4).*t));
end


subplot(4,1,1),plot(time,xsp(1,:));ylabel('\^u');title('Phugoid Mode');
subplot(4,1,2),plot(time,xsp(2,:));ylabel('\alpha');
subplot(4,1,3),plot(time,xsp(3,:));ylabel('\^q');
subplot(4,1,4),plot(time,xsp(4,:));ylabel('\theta');xlabel('time');

%trajectory
figure;
sigma = real(eval(3,3));
omega = imag(eval(3,3));

period = abs(2*pi/omega);
time = linspace(0,6*period,201);
vec = evec;
for i = 1:length(time)
   t = time(i);
   factor = 2*exp(sigma*t)/(sigma^2+omega^2);
   u(i) = 2*exp(sigma*t)*(real(vec(1,3))*cos(omega*t)-imag(vec(1,3))*sin(omega*t));
   w(i) = 2*exp(sigma*t)*(real(vec(2,3))*cos(omega*t)-imag(vec(2,3))*sin(omega*t));
   theta(i) = 2*exp(sigma*t)*(real(vec(4,3))*cos(omega*t)-imag(vec(4,3))*sin(omega*t));
   xe(i) = u0*t + factor*(real(vec(1,3))*(sigma*sin(omega*t)-omega*cos(omega*t))+...
      imag(vec(1,3))*(sigma*cos(omega*t)+omega*sin(omega*t)));
   ze(i) = factor*((real(vec(2,3))-u0*real(vec(4,3)))*(sigma*sin(omega*t)-omega*cos(omega*t))+...
      (imag(vec(2,3))-u0*imag(vec(4,3)))*(sigma*cos(omega*t)+omega*sin(omega*t)));
end

subplot(2,1,1);plot(xe,ze);
title('Aircraft Trajectory.  Phugoid Mode.')
ylabel('z_E');
subplot(2,1,2);plot(xe,u);title('Disturbance Velocity');ylabel('u');xlabel('x_E');
figure
plot(xe-u0*time,ze);axis equal;title('Trajectory in Moving Coordinate System');