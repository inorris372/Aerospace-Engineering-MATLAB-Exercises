%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Test 2 AerE 331          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1
%Part A

% denom1 = [1 2.21 1.78 2.13 0.07];
% roots(denom1)
% 
% %Part B
% k = 2.77;
% num1 = [.73 1.43 .3 .24];
% trans1 = tf(num1,denom1);
% trans1cl = feedback(k*trans1,1);
% a = 0:.001:100;
% figure(1)
% rlocus(trans1)
% grid
% figure(2)
% step(trans1cl,a)
% grid
% %% Problem 2
% %Part B
% 
% 
% trans2 = tf(40000, [1 51 1050 1000 40000]);
% CL = feedback(trans2, 1);
% step(CL) % Step response
% t = 0:.01:5;
% figure(3)
% lsim(CL,t,t) % Ramp response
% 
% %Part D
% 
% 
% trans2 = tf(40000, [1 51 1050 1000 40000]);
% CL = feedback(trans2, 1);
% step(CL) % Step response
% t = 0:.01:5;
% figure(4)
% lsim(CL,t,t) % Ramp response





%% Problem 3
%Part A

num2 = 20000;
% denom2 = [3 153 3150 3000];
% trans3 = tf(num2,denom2);
% [Gm,Pm,Wgm,Wpm] = margin(trans3);
% roots(denom2)
% 
% figure(5)
% bode(trans3)
% margin(trans3)
% grid

%Part B
e = 18.66025;
d = 1.339746;
numl = [1 d]*e;
denoml = [1 e]*d;
k = 10^(-13.7/20);
lag = tf([1000 1],[1000*k 1]);
lead = tf(numl,denoml);
num2a = num2;
denom2a = [3 153 3150 3000];
Gopen = tf(num2a,denom2a);
Gaug1 = Gopen*lead;
Gaug2 = Gopen*lead*lag;
figure(55)
bode(Gaug1)
grid
figure(56)
Y = feedback(Gaug2,1);
step(Y)
t2 = 0:.01:5;
grid
figure(57)
lsim(Y,t2,t2)
grid



% %% Problem 4
% % PROBLEM 4 "The Importance of FRF Phase"
% % 4(a):
% %Plant: Gp =K/s*(tau*s+1);
% %==============================
% %time delay each way:
% dmin = 0.02042;                       % ***Minimum Time Delay***
% dmax_sf = 0;                          % Maximum Time Delay Scale Factor
% d = dmin*(1 + dmax_sf*rand(1,1))      % Actual Random Delay (sec)
% tau=.5;                               % Drone Time Constant
% K = 100;                              % Open Loop Static Gain
% G = tf(K,[tau 1 0],'iodelay',d);      % Forward Loop TF w/o Compensation
% H = tf(1,1,'iodelay',d);              % Feedback Loop TF
% % UNITY STATIC GAIN LEAD COMPENSATION:
% wgc = 14;                             % wgc prior to Lead
% PHmax = 50;                           % Chosen Compensator Max. Phase (Deg)
% phmax = PHmax*pi/180;                 % Convert PHmax to radians
% rw = (1+sin(phmax))/(1-sin(phmax));   % Computed to = -17.5dB
% wc = 23.4;                            % Place phimax @ -8.75dB = -(rw dB)/2
% w1 = wc/sqrt(rw);                     % Compensator Low Freq. (rad)
% w2 = rw*w1;                           % Compensator Hight Freq. (rad)
% Gc = tf([(1/w1) 1],[(1/w2) 1]);       % Compensator Transfer Function
% GG = Gc*G;                            % Updated Forward Loop TF
% GGH = GG*H;                           % Updated Open Loop TF
% figure(12)
% w=logspace(-1,3,1000);
% bode(GGH,w)
% [Gm,Pm,Wgm,Wpm] = margin(GGH)
% margin(GGH)
% title('Open Loop Bode Plot for K=100 & Lead Compensation')
% grid
% 
% WW = feedback(GG,H);                  % Closed Loop TF
% figure(13)
% step(WW)
% grid
% title('Unit Step Response for K=100 & Lead Compensation')
% 
% %PROBLEM 4: “The Influence of Random Phase Distortion”
% Fs = 44100;                            % Digital Sampling Frequency (Hz)
% load ch1                               % Music file (ch1.mat)
% sound(ch1,Fs)                          % Listen to undistorted music clip
% pause
% CH1 = fft(ch1);                        %Fourier Transform of ch1
% N = length(ch1);
% df = Fs/N;
% f = 0:N-1;
% f = df*f';                             %Frequencies associated with CH1
% CH1PH = angle(CH1);                    %Phase associated with CH1
% figure(21)
% plot(f,(180/pi)*CH1PH)
% title('FFT Phase of CH1')
% grid
% hold on
% %=================================
% % ADD RANDOM PHASE DISTORTION
% N2 = fix(N/2);
% sigma_dPH = 1; % *** SELECTABLE std. deviation of random distortion ***
% dPH = sigma_dPH*randn(N2,1);            %Random Phase Distortion
% %dPH = (pi/2)*sin(2000*(2*pi)*f2/f2max);
% f2 = f(1:N2);
% f2max = f(N2);
% dPH = [0 ; dPH ; -flipud(dPH) ];
% CH1PHH = CH1PH + dPH;                  %New Phase of CH1
% plot(f,(180/pi)*CH1PHH,'r')
% title('Comparison of Original & Modified Phase of CH1')
% hold off
% pause
% CH11 = abs(CH1).*exp(1i*CH1PHH);       %New CH1   
% ch11 = real(ifft(CH11));               %New ch1
% sound(ch11,Fs)
% figure(22)
% plot(ch1)
% hold on
% plot(ch11,'r')
% title('Original & Phase Distorted Music')
% grid
% hold off

