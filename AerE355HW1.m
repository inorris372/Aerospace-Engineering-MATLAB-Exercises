%%%%%%%%%%%%%%%%%%%%%%%%%
%  AERE 355 HW1         %
%  Ian Norris           %
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%  Function of Range(Cl) when Cl varies between 0.1 and 1.1 %%%

Cl = [.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.1]; %Coefficient of Lift
Cd0 = .035;  %Zero Drag Coefficient
e = .85;  %Induced Drag
p = .07651;  %Air Density
W = 200;  %Weight of Aircraft
AR = 8;  %Aspect Ratio
Cd = e + (Cl.*Cl)./(pi * e * AR);  %Expression for Coefficient of Drag
alpha = atan(-Cd./Cl);  %Expression for Alpha - Angle of Attack
v = sqrt((2 * W)./(AR * p .* cos(alpha)- Cd .* p .* sin(alpha)));
     % Expression for Airspeed Velocity of the Aircraft
vx = sqrt((v.*v)./((tan(alpha)).^2 + 1));  %Horizontal Velocity Vector
range = -2000./tan(alpha);  %Range of Aircraft - how far it can go

figure;
plot(Cl,range)  %Graphs Coefficient of Lift vs. Range of Aircraft

figure(2);
plot(Cl,vx)  %Graphs Coefficient of Lift vs. Velocity of Aircraft
