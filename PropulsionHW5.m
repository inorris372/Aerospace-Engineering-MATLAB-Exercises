%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Ian Norris                     %
%          Propulsion Assignment # 5        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear home
%% Problem 2
%% Part A
gamma = 1.4;                         %[Unitless] Ratio of Specific Heats
T0 = 220;                            %[K] Initial Temperature
cp = 1005;                           %[Unitless] Coefficient of Pressure
hf = 4.42e7;                         %[J/kg] Heating value of Fuel
M0 = 2.5;                            %[Unitless] Initial Mach Number
tlam = 7;                            %[Unitless] Tau Lambda
tlamAB = linspace(7,9,100);          %[Unitless] Tau Lambda for Afterburner

a0 = sqrt((gamma-1)*cp*T0);          %[m/s] Velocity of Sound
tr = 1 + ((gamma-1)/2)*M0^2;         %[Unitless] Tau Tt0/T0
ST = a0.*(sqrt((2.*tlamAB./(gamma-1).*(1-4*tlam./(tlam+tr).^2)))-M0);
%[N/(kg/s)] Specific Thrust
f = cp.*T0./hf.*(tlamAB-tr);         %[Unitless] Fuel to Air Ratio
TSFC = f./ST;                        
%[(kg/s)/N] Thrust Specific Fuel Consumption

figure(1)
plot(tlamAB,ST)
xlabel('Tau Lambda AB')
ylabel('Specifc Thrust [N/(kg/s)]')

figure(2)
plot(tlamAB,TSFC)
xlabel('Tau Lambda AB')
ylabel('TSFC [(kg/s)/N]')

%% Part B

tlamAB = 8;                          %[Unitless] Tau Lambda for Afterburner
tlam = linspace(6.5,7.5,100);        %[Unitless] Tau Lambda

ST = a0.*(sqrt((2.*tlamAB./(gamma-1).*(1-4*tlam./(tlam+tr).^2)))-M0);
%[N/(kg/s)] Specific Thrust
f = cp.*T0./hf.*(tlamAB-tr);         %[Unitless] Fuel to Air Ratio
TSFC = f./ST;          
%[(kg/s)/N] Thrust Specific Fuel Consumption

figure(3)
plot(tlam,ST)
xlabel('Tau Lambda AB')
ylabel('Specific Thrust [N/(kg/s)]')

figure(4)
plot(tlam,TSFC)
xlabel('Tau Lambda AB')
ylabel('TSFC [(kg/s)/N]')

%% Problem 3
%% Part A

M0 = linspace(0,3.5,100);            %[Unitless] Initial Mach Number
gamma = 1.4;                         %[Unitless] Ratio of Specific Heats
pic = 10;                            %[Unitless] Pi of Compressor
T0 = 210;                            %[Kelvin] Initial Temperature
cp = 1005;                           %[Unitless] Coefficient of Pressure
hf = 4.42e7;                         %[J/kg] Heating value of Fuel
tlam = 7;                            %[Unitless] Tau Lambda
tlamAB = 8;                          %[Unitless] Tau Lambda for Afterburner
tc = (pic)^((gamma-1)/gamma);        %[Unitless] Tau Compressor
tr = 1 + ((gamma-1)/2).*M0.^2;       %[Unitless] Tau Tt0/T0
tt = 1 - tr./tlam.*(tc-1);           %[Unitless] Tau Turbine

ST = a0.*(sqrt(2.*tlamAB./((gamma-1).*tt.*tc.*tr).*(tt.*tc.*tr-1))-M0);
%[N/(kg/s)] Specific Thrust
f = cp.*T0./hf.*(tlamAB-tr);         %[Unitless] Fuel to Air Ratio
TSFC = f./ST;
%[(kg/s)/N] Thrust Specific Fuel Consumption

figure(5)
plot(M0,ST)
xlabel('Initial Mach Number')
ylabel('Specific Thrust [N/(kg/s)]')

figure(6)
plot(M0,TSFC)
xlabel('Initial Mach Number')
ylabel('TSFC [(kg/s)/N]')

%% Part B

M0 = linspace(0,5,1000);             %[Unitless] Initial Mach Number
tr = 1 + ((gamma-1)/2).*M0.^2;       %[Unitless] Tau Tt0/T0
tt = 1 - tr./tlam.*(tc-1);           %[Unitless] Tau Turbine
ST = a0*(sqrt(2*tlamAB./((gamma-1).*tt.*tc.*tr).*(tt.*tc.*tr-1))-M0);
%[N/(kg/s)] Specific Thrust

for i = 1:999
    mini = abs(ST(i+1));
    if(mini<ST(i))
        keep = i+1;
    end
end

sprintf('The maximum Mach number at which these engines will be able to\n')
sprintf('fly is %d.',keep*5/1000)
%% Problem 4
%% Part A

v = 900;                          %[m/s] Velocity
h = 10;                           %[km] Height
d = 0.41351;                      %[kg/m^3] Air Density
st = 223.26;                      %[K] Air Static Temperature
hf = 4.4e7;                       %[J/kg] Heating Value for the Fuel
nta = 0.5;                        %[m^2] Nozzle Throat Area
stagt = 1500;                     %[K] Stagnation Temperature Max.
R = 287;                          %[J/(kg*K)] Gas Constant
gamma = 1.4;                      %[Unitless] Ratio of Specific Heats
cp = 1005;                        %[Unitless] Pressure constant
T0 = st;                          %[K] Initial Temperature
Tt4 = stagt;                      %[K] Temperature downstream of burner

a0 = sqrt(gamma*R*st);            %[m/s] Speed of Sound
Mach = v/a0;                      %[Unitless] Mach Number
tlamb = Tt4/T0;                   %[Unitless] Tau lambda Tt4/T0

sprintf('M0 = %g and Tau Lambda = %g',Mach,tlamb) 

%% Part B

tr = 1 + (gamma-1)/2*Mach^2;      %[Unitless] Tt0/T0
ST = a0*Mach*(sqrt(tlamb/tr)-1);  %[N/(kg/s)] Specific Thrust
P0 = R*d*T0;                      %[N/m^3] Initial Pressure
Pt8 = tr^(gamma/(gamma-1))*P0;    %[N/m^3] Pressure after turbine
Gam = sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)));
mdot = Pt8*Gam*nta/sqrt(R*(stagt+15.09));    %[kg/s] Mass flow
Thrust = ST*mdot;                            %[N] Thrust
sprintf('Thrust = %g N', Thrust)

%% Part C
f = cp*T0/hf*(tlamb-tr);          %[Unitless] Fuel to Air Ratio
mdotf = mdot*f;                   %[kg/s] fuel mass flow rate
sprintf('mdotf = %g kg/s',mdotf)







































