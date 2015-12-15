%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AerE 331 Hw-6 Ian Norris      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Problem 2
k = 1:1:10;
hz = 2.^k;
% figure(1)
% plot(k,hz)
% axis([0 10 0 1200])

%% Problem 3

%Part E

w = .01:.01:100;
T = .1;
Hs = (1 + w.*1i).^(-1);
HsdB = 20*log10(abs(Hs));
HsP = (180/pi)*angle(Hs);
z1 = exp(-1i*w*T);

figure(2)
semilogx(w,HsdB)
hold on
plot(w,z1)
xlabel('Frequency')
ylabel('dB')
hold off

figure(7)
plot(w,HsP)
xlabel('frequency')
ylabel('phase')

%% Problem 4
f = 20;
t = 0:.05*10^-3:.05;
x = sin(2*pi*f*t);
t2 = 0:.005:.05;
x2 = sin(2*pi*f*t2);
figure(3)
plot(t,x)
hold on
plot(t2,x2)
xlabel('time (seconds)')
ylabel('frequency')
hold off

%% Problem 5

%Part B
T = 0.2;
K=5;
s = tf('s');
z1 = exp(-T*s);
Gcz = K*(1-z1)/T;
Ghats = Gcz*(1/T)*(1-z1)/s;
Gcs =K*s;
bode(Gcs,'b', Gcz,'r')
step(Ghats)
title('Analog (blue) & digital (red) D-controller Bode plots')
grid
