%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Ian Norris                      %
%        AerE 411 Propulsion             %
%        Assignment #4                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
%% Problem 1
%1. (1 point) (a) Calculate and carefully plot the specific thrust, F/m? 
% (in N.s/kg), and thrust specific fuel consumption, TSFC (in kg/N.s), vs 
% design pressure ratio, ?c for an ideal turbojet. The flight Mach number, 
% M0 = 2.7, and you are given ? = 1.4, T0 = 400 K, hf = 4.4194 × 107 J/kg, 
% cp = 1004.5 J/kg.K, and ?? = 6.8.  Plot the results over the range of 
% ?c = 1 to where the thrust goes to zero.  (b) Obtain an analytic 
% expression for ?c that leads to zero thrust, in terms of ?, M0, and ??. 
% Check that the value calculated using the data of part (a) agrees with 
% the graphical result of part (a).

%Part  A
T0 = 400;                      %[Kelvin]  Initial Temperature 
tlamb = 6.8;                   %[Unitless]
gamma = 1.4;                   %[Unitless]
M0 = 2.7;                      %[Unitless] Initial Mach Number
cp = 1004.5;                   %[J/kg.K] Coefficient of Pressure
R = 287;                 %[Nm/molK] Ideal Gas Constant
hf = 4.4194e7;                 %[J/kg] Final Specific Enthalpy
a0 = sqrt(gamma*R*T0);         %[m/s] Initial Speed of Sound
pic = linspace(1,40,1000);     %[Unitless]  Compressor Pressure Ratio
tr = 1 + (gamma-1)/2*M0^2;     %[Unitless]  
tc = pic.^((gamma-1)/gamma);   %[Unitless]  Compressor Temp. Ratio
tt = 1-tr./tlamb.*(tc-1);      %[Unitless]Temperature Ratio (Final/Initial)
ST = a0.*(sqrt(2./(gamma-1).*tlamb./(tc.*tr).*(tt.*tc.*tr-1))-M0);  
%[F/(kg/s)]  Specific Thrust 

%Part B
f = cp.*T0./hf.*(tlamb-tc.*tr);
TSFC = f./ST;

figure(1)
plot(pic,TSFC*10e6,'r')
hold on
plot(pic,ST,'b')
xlabel('Compressor Ratio [Unitless]')
ylabel('Specific Thrust,TSFC - scaled(10e6)')
axis([1 40 0 1000])
legend('TSFC','ST')

%% Problem 2
%2. (1 point) Calculate the specific thrust, F/m? and thrust specific fuel 
% consumption, TSFC for a ramjet with the conditions as given for the 
% turbojet in Problem 1, except that ?? = 8.0.

tlamb = 8;                     %[Unitless]
tt = 1-tr./tlamb.*(tc-1);      %[Unitless]Temperature Ratio (Final/Initial)
T9 = T0.*tlamb./(tc.*tr);      %[Kelvin]  Final Temperature
ST = a0.*(sqrt(2./(gamma-1).*tlamb./(tc.*tr).*(tt.*tc.*tr-1))-M0);  
%[F/(kg/s)]  Specific Thrust 

%Part B
f = cp.*T0./hf.*(tlamb-tc.*tr);   %Fuel Consumption
TSFC = f./ST;
display(TSFC(1000))
%Thrust Specific Fuel Consumption

figure(2)
plot(pic,TSFC*10e6,'r')
hold on
plot(pic,ST,'b')
xlabel('Compressor Ratio [Unitless]')
ylabel('Specific Thrust,TSFC - scaled(10e6)')
axis([1 40 0 1000])
legend('TSFC','ST')

%% Problem 3
%3. (1 point) Design a nozzle for an ideal rocket that has to operate at 
% 25 km altitude (pa = 2.549 k Pa) and give 5000 N thrust at a chamber 
% pressure of pc = 2.068 M Pa and a chamber temperature of Tc = 2800 K. 
% Assuming gamma = 1.30, and R = 355.4 J/kg/K, determine the throat area, 
% exit area, throat velocity, and exit temperature.

h = 25;
Pa = 2.549*1000;
thrust = 5000;
Pc = 2.068*10^6;
Tc = 2800;
roc = Pc/(R*Tc);
Tstar = Tc;
gamma = 1.3;
R = 355.4;
Mt = 1;
at = sqrt(Mt*gamma*R);



TstarbyTe = 1 + (gamma-1)/2*Mt^2;
Te = Tstar/TstarbyTe;      %[Kelvin]  Final Temperature
ae = sqrt(gamma*R*Te);
Pe = 1/((TstarbyTe)^(1/(gamma-1))/Pc);
roe = Pe/(R*Te);
ue = sqrt(2*gamma*R/(gamma-1)*Tc*(1-(Pe/Pc)^((gamma-1)/gamma)));
Me = ue/ae;
AstarbyAe = ((gamma+1)/2)^(1/(gamma-1))*(Pe/Pc)^(1/gamma)*sqrt((gamma+1)...
    /(gamma-1)*(1-(Pe/Pc)^((gamma-1)/gamma)));
Cf = sqrt(2*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1))*(1-...
    (Pe/Pc)^((gamma-1)/gamma)))+(Pe/Pc-Pa/Pc)*1/AstarbyAe;
Astar = thrust/(Pc*Cf);
Ae = Astar/AstarbyAe;
mdote = (thrust-(Pe-Pa)*Ae)/ue;

display(ue)
display(Te)
display(Astar)
display(Ae)


%% Problem 4
%4. (1 point) Consider a rocket with a translatable skirt, such that it has
% two design altitudes, 10,000 m and 20,000 m.  The ambient pressure versus
% altitude may be approximated by the formula, p0 = pSL exp(-h/7000), where 
% pSL = 1.01×105 Pa is the sea-level standard pressure. For this rocket the 
% Summerfield criterion is ps = pa/2.718.  Take gamma = 1.22 and the chamber 
% pressure, pc = 50 × pSL . The skirt is translated at just the altitude 
% where the CF for the skirt in the 10,000 m design condition equals the CF
% for the skirt in the 20,000 m design condition.  Calculate and carefully 
% plot CF over the range 0 ? h ? 40, 000 m. Carefully locate the altitudes, 
% htrans, hD1 and hD2. Include on the graph the envelope of a perfectly 
% expanded rocket.

TstarbyTe = 1 + (gamma-1)/2*Mt^2;

h = linspace(0,40000,80000);
Dalt1 = 10000;
Dalt2 = 20000;
pSL = 1.01e5;
Pa = pSL.*exp(-h./7000);
Pc = 50*pSL;
%Pe = 1/((TstarbyTe)^(1/(gamma-1))/Pc);
Pe = Pa;
Pe2 = Pa(Dalt1*2);
Pe3 = Pa(Dalt2*2);
gamma = 1.22;

AstarbyAe = ((gamma+1)/2)^(1/(gamma-1)).*(Pe./Pc).^(1/gamma).*...
    sqrt((gamma+1)/(gamma-1)*(1-(Pe./Pc).^((gamma-1)/gamma)));
Cf = sqrt(2.*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1)).*(1-...
    (Pe./Pc).^((gamma-1)/gamma)))+(Pe./Pc-Pa./Pc).*1./AstarbyAe;

AstarbyAe2 = ((gamma+1)/2)^(1/(gamma-1)).*(Pe2./Pc).^(1/gamma).*...
    sqrt((gamma+1)/(gamma-1)*(1-(Pe2./Pc).^((gamma-1)/gamma)));
Cf2 = sqrt(2.*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1)).*(1-...
    (Pe2./Pc).^((gamma-1)/gamma)))+(Pe2./Pc-Pa./Pc).*1./AstarbyAe2;

AstarbyAe3 = ((gamma+1)/2)^(1/(gamma-1)).*(Pe3./Pc).^(1/gamma).*...
    sqrt((gamma+1)/(gamma-1)*(1-(Pe3./Pc).^((gamma-1)/gamma)));
Cf3 = sqrt(2.*gamma^2/(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1)).*(1-...
    (Pe3./Pc).^((gamma-1)/gamma)))+(Pe3./Pc-Pa./Pc).*1./AstarbyAe3;

figure(3)
plot(h,Cf,'b')
hold on
plot(h,Cf2,'r')
hold on
plot(h,Cf3,'g')
hold off
xlabel('Height(Meters)')
ylabel('Cf')
axis([0 40000 1.4 2.1])
legend('Cf','Cf2','Cf3')


