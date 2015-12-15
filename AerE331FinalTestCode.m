%%   Ian Norris - AerE 331 Final Exam

%% Problem 1 

%Part C
Gc = 1;
Gs = tf([10],[1 10]);
Gp = tf([50],[1 1.4 0]);
Gopen = Gc*Gs*Gp;
num1 = 500;
denom1 = [1 11.4 14 0];
Gclosed = Gopen/(Gopen+1);
[A,B,C,D] = tf2ss(num1,denom1);

Poles = [-2 -2+1.5i -2-1.5i];
ks = place(A,B,Poles)

%Part D
figure(1)   %Closed loop step response for ks = 1
step(Gclosed)

figure(2)
Gopen2 = Gopen*.025;
Gclosed2 = feedback(Gopen2,tf(1,1));
step(Gclosed2)

%% Problem 2

%Part A
A = [-12 -400; 1 0];
lambda1 = eig(A)

%Part D
om = 0:.01:100;
w = 0.1;
u = (om.*(1i)).^2+12*(om.*(1i))+400;
h = w./u;
sw = 10*log10(abs(h).^2);

figure(3)
plot(om,sw)
xlabel('frequency (w)')
ylabel('dB')

%Part E
numC = 0.1;
denomC = [1 12 400];
transC = tf(numC,denomC);
r = c2d(transC,2*pi/1000)

u = rand(1,5000);
x = linspace(0,100,length(u));
y(1)=0;
y(2)=1.923*1e-6*u(1);

for k = 3:5000
    y(k) = 1.912*y(k-1)-0.9274*y(k-2)+1.923*1e-6*u(k-1)+1.875*1e-6*u(k-2);
end

figure(77)
title('Simulated Turbulence')
plot(x,y)



%Part F
figure(4)
psd(sw)
%% Problem 3

%Part B
num3 = [2 0 0 0];
denom3 = [1 -2.5 2.26 -0.74];
trans3 = tf(num3,denom3,.05);
roots(denom3)  


%Part D
trans3d = feedback(trans3,1);
figure(5)
rlocus(trans3d)

%Part E
figure(6)
impulse(trans3,'b')
hold on
impulse(trans3d,'g')

%% Problem 4

%Part B
s = tf('s');
%e = 2.7182818;
z = exp(-s.*.02);
ZOH = (1-z)/s;
Gz = .004/(1-.98*z);
Gp = tf(4.6,[1 0.76 4.55]);
Gs = tf(.2,[1 1]);
Gc = Gp*Gz*ZOH*(1/.02)

wA = feedback(Gs*Gp,1);
wD = feedback(Gs*Gc,1);

%Part C
figure(76)
step(wA,'b')
hold on
step(wD,'r')
grid
hold off

%Part D
figure(75)
bode(wA,'b')
hold on
bode(wD,'r')
grid
hold off


% %% Problem 5
% 
% %PROGRAM NAME: ButterAveMaria.m
% %This code carries out a decomposition of the vocal portion of the audio
% % track (ch1.mat MUST be present in the Matlab workspace.)
% %======================================================================
% % Analog & Digital Filter Design & Properties:
% fs = 44200;     % Sample frequency in Hz
% Ts = 1/fs;      % Sampling period
% fN = fs/2;      % Nyquist Frequency
% f12=[400 4000]; %  Vocals BP Filter (Hz)
% F12 = f12/fN;   % Normalized BW
% Norder = 4;     % Because the filter is BP, its order will = 2*Norder
% [Bs,As] = butter(Norder , F12,'s'); %Max frequency possible is 1 RAD/SEC
% %ANALOG FILTER:
% Gs=tf(Bs,As);
% s=tf('s');
% npwr=length(As);
% BBs=0; AAs=0;
% % Convert frequencies to REAL frequencies
% c=(2*pi*fN)^-1;   % wnorm = c*wreal
% for k=1:npwr
%     BBs = BBs + Bs(k)*(c*s)^(npwr-k);
%     AAs = AAs + As(k)*(c*s)^(npwr-k);
% end
% GGs = BBs/AAs; % Analog filter re: REAL frequencies
% % DIGITAL FILTER:
% [B,A] = butter(Norder,F12);
% Gz = tf(B,A,Ts);
% %Comparison of Analog & Digital FRS & Impulse Responses:
% figure(1)
% bode(GGs,'k',Gz,'r')
% grid
% title('Analog (black) & Digital (red) BP Filter Bode Plots')
% pause   %PAUSE 1
% figure(2)
% impulse(GGs,'k',Gz,'r')
% grid
% title('Analog (black) & Digital (red) BP Filter Impulse Responses')
% pause   %PAUSE 2
% %=================================================================
% % Carry out the Vocals Filtering:
% %-----------------------------------------------
% % ***** Run ch1 through the digital filter *****
% nch1 = length(ch1);
% ch1F = zeros(nch1,1);
% tvec=Ts:Ts:nch1*Ts;
% nA=length(A); nB=length(B);
% for k = nA+1:nch1
%     ch1F(k)=-A(2:nA)*ch1F(k-1:-1:k-nA+1) + B*ch1(k:-1:k-nB+1);
% end
% figure(3)  %+++ ZOOM IN TO COMPARE UNFILTERED & FILTERED IN DETAIL +++
% plot(tvec,ch1)
% grid
% title('Unfiltered Track')
% xlabel('Time (sec)')
% sound(ch1,fs)
% pause   %PAUSE 3
% hold on
% plot(tvec,ch1F,'r')
% title('Unfiltered (blue) & Once-Filtered (red) Tracks')
% sound(ch1F,fs)  % Play filtered ch1
% pause % PAUSE 4
% chdiff=ch1-ch1F;
% figure(4)
% plot(tvec,chdiff)
% title('Track Residual After Once-Filtered Vocals Have been Subtracted')
% xlabel('Time (sec)')
% grid
% sound(chdiff,fs)  % Play Residual after removing ch1F from the track
% pause   % PAUSE 5
% %---------------------------------------------
% % ***** Run Ch1F backwards through the filter *****
% ch1Fi = flipud(ch1F);
% ch1F2i = zeros(nch1,1);
% for k = nA+1:nch1
%     ch1F2i(k)=-A(2:nA)*ch1F2i(k-1:-1:k-nA+1) + B*ch1Fi(k:-1:k-nB+1);
%     
% end
% ch1F2 = flipud(ch1F2i);
% figure(5)
% plot(tvec,ch1)
% hold on
% plot(tvec,ch1F2,'r') %++ZOOM IN TO COMPARE UNFILTERED & FILTERED IN DETAIL++
% title('Unfiltered (blue) & Twice-Filtered (red) Tracks')
% xlabel('Time (sec)')
% grid
% chdiff2 = ch1 - ch1F2;
% sound(chdiff2,fs)
% figure(6)
% plot(tvec,chdiff2)
% grid
% title('Track Residual After Twice-Filtered Vocals Have been Subtracted')
% xlabel('Time (sec)')






















