%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Ian Norris                %
%       Astrodynamics II          %
%       Problem Set 5             %
%       Gravity-Turn Rocket       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
global m0 g0 T A Cd rh0 H0 Re hgr_turn tf md
% Launch Site: Pad A at Kennedy Space Center
Alt = 14.6;              %[m] Alt above sea level
m_star = 3100;           %[kg] Payload Mass
m_prop1 = 88365;      % [kg] Propellant mass
m_prop2 = 106261 + 629340;  % [kg] Propellant mass
eps1 = .11;           % Stage 1 Structure factor 
eps2 = .13;           % Stage 2 Structure factor
m_st1 = eps1*m_prop1/(1-eps1);
m_st2 = eps2*m_prop2/(1-eps2);
m2f = m_prop2 + m_star + m_st2;
m1f = m_prop1 + m2f + m_st1;

mu = 3.986e5;  % Gravitational Parameter
g0     = 9.81;        % [m/s^2] Constant at its sea-level value
Isp1   = 280 ;        % [s]  Specific Impulse Vega P80 Rocket
Isp2   = 450 ;        % [s]  Spefific Impulse Space Shuttle Main Engines
% ac = -Beta*Ve/m-g; %Initial Acceleration
rp = 35780 + 6378;    % [km] Radius to Perigee 
ra = 72 + 6378;       % [km] Radius to Apogee
a = (rp+ra)/2;        % [km] Semi-Major Axis
ee = ra/a-1;          % Eccentricity of Goal Orbit
vbo = 15*pi/180;      % [radians] true anomaly angle at burn out
rentrance = a*(1-ee^2)/(1+ee*cos(vbo)); % [km] Radius @ GTO Insert
Vb = sqrt(mu*(2/rentrance-1/a));  %Velocity @ GTO Insert
fa = 180*atan(ee*sin(vbo)/(1+ee*cos(vbo)))/pi; %Flight Path Angle @ GTO Insert
V0 = 0;
Vesh = Isp2*g0;          % [m/s] Space Shuttle Exhaust velocity
Vveg = Isp1*g0;          % [m/s] Vega P80 Exhaust Velocity
delVa = Vveg*log(m1f/(m1f-m_prop1)); % [m/s] Delta V induced by 1st stage
delVb = Vesh*log(m2f/(m2f-m_prop2)); % [m/s] Delta V induced by 2nd stage
delV = delVa + delVb;
m0  = m1f;            % [kg] Initial mass
A   = 20;             % [m^2]Frontal area
Cd  = 0.14 ;          % Drag coefficient = constant value
rh0 = 1.2;            % [kg/m^3] Air Density at Sea-Level
H0 = 8000;            % [m] Density scale height
Re = 6378e3;          % [m] Earth's radius
hgr_turn = 200;       % [m] Rocket starts gravity turn when h = hgr_turn

mf = m0 - m_prop1;  % [kg] Final mass of the rocket
t0 = 0;             % Rocket launch time
% Launch initial conditions:
gamma0 = 89.5/180*pi;       % Initial flight path angle
v0 = 0;   % Velocity (m/s)  % Earth's Rotation considered in eq of motion.
x0 = 0;   % Downrange distance [km]
h0 = Alt; % Launch site altitude [km]
vD0 = 0;  % Loss due to drag (Velocity)[m/s]
vG0 = 0;  % Loss due to gravity (Velocity)[m/s]
state0   = [v0, gamma0, x0, h0, vD0, vG0];
beta1 = -21400;
beta2 = -120;
         

%tburn1 = 106.8;      % [s] Fuel burn time, first stage
tburn1 = -(m1f-m_star)*(1-eps1)/beta1; % [s] Fuel burn time, first stage
%tburn2 = 480;        % [s] Fuel burn time, second stage
tburn2 = -(m2f-m_star)*(1-eps2)/beta2; % [s] Fuel burn time, second stage
md1 = (m_prop1)/tburn1;  % [kg/s]Propellant mass flow rate (1st stage)
md2 = (m_prop2)/tburn2;  % [kg/s]Propellant mass flow rate (2nd stage)
T1   = md1*(Vesh);    % [N] Thrust (mean)
T2   = md2*(Vveg);    % [N] Thrust (mean)

for g = 1:2        
    if g == 1
        T = T1;
        md = md1;
        m0 = m1f;
        tf = t0 + tburn1;  % The time when propellant is completely 
        %burned and the thrust goes to zero
        t_range = [t0,tf];  % Integration interval
        % Solve initial value problem for ordinary differential equations
        [t,state] = ode45(@RocketDynEq,t_range,state0) ;
        v1  = state(:,1)/1000;     % Velocity [km/s] 
        gamma1 = state(:,2)*180/pi;    % Flight path angle  [deg]
        x1  = state(:,3)/1000;     % Downrange distance [km]
        h1  = state(:,4)/1000;     % Altitude[km]
        vD1 = -state(:,5)/1000;    % Loss due to drag (Velocity)[m/s]
        vG1 = -state(:,6)/1000;    % Loss due to gravity (Velocity)[m/s]
     
        % VEGA Rocket: First Stage P80

        fprintf('\n VEGA Rocket: First Stage P80\n')
        fprintf('\n Propellant mass           = %4.2f [kg]',m_prop1)
        fprintf('\n Gross mass                = %4.2f [kg]',m1f)
        fprintf('\n Isp                       = %4.2f [s]',Isp1)
        fprintf('\n Thrust(mean)              = %4.2f [kN]',T/1000)
        fprintf('\n Initial flight path angle = %4.2f [deg]',gamma0*180/pi)
        fprintf('\n Final speed               = %4.2f [km/s]',v1(end))
        fprintf('\n Final flight path angle   = %4.2f [deg]',gamma1(end))
        fprintf('\n Altitude                  = %4.2f [km]',h1(end))
        fprintf('\n Downrange distance        = %4.2f [km]',x1(end))
        fprintf('\n Drag loss                 = %4.2f [km/s]',vD1(end))
        fprintf('\n Gravity loss              = %4.2f [km/s]',vG1(end))
        fprintf('\n');
            
        figure(1)
        plot(t,h1,'b');
        hold on;
        grid on;
        %plot(t,h,'.b');
        xlabel('time[s]');
        ylabel('Altitude[km]');

        figure(2)
        plot(t,x1,'b');
        hold on;
        grid on;
        %plot(t,x,'.b');
        xlabel('time[s]');
        ylabel('Downrange Distance[km]');
    
        figure(3)
        plot(t,gamma1,'b');
        hold on;
        grid on;
        %plot(t,gamma,'.b');
        xlabel('time[s]');
        ylabel('gamma [radians]');
        
    else
        hgr_turn = h1(end);
        m0 = m2f;
        T = T2;
        md = md2;
        t0 = tf;
        tfb = tf + tburn2;
        t_range2 = [t0,tfb];
        % Solve initial value problem for ordinary differential equations
        state2 = [v1(end), gamma1(end), x1(end), h1(end), vD1(end), vG1(end)];
%         [t,state] = ode45(@RocketDynEq,t_range2,state2) ;
%         v2  = state(:,1)/1000;     % Velocity [km/s] 
%         gamma2 = state(:,2)*180/pi;    % Flight path angle  [deg]
%         x2  = state(:,3)/1000;     % Downrange distance [km]
%         h2  = state(:,4)/1000;     % Altitude[km]
%         vD2 = -state(:,5)/1000;    % Loss due to drag (Velocity)[m/s]
%         vG2 = -state(:,6)/1000;    % Loss due to gravity (Velocity)[m/s]
        
        % Launch System using Space Shuttle Main Engine
        fprintf('\n Launch System using Space Shuttle Main Engine\n')
        fprintf('\n Propellant mass           = %4.2f [kg]',m_prop2)
        fprintf('\n Gross mass                = %4.2f [kg]',m2f)
        fprintf('\n Isp                       = %4.2f [s]',Isp2)
        fprintf('\n Thrust(mean)              = %4.2f [kN]',T/1000)
        fprintf('\n Initial flight path angle = %4.2f [deg]',gamma0*180/pi)
%         fprintf('\n Final speed               = %4.2f [km/s]',v2(end))
%         fprintf('\n Final flight path angle   = %4.2f [deg]',gamma2(end))
%         fprintf('\n Altitude                  = %4.2f [km]',h2(end))
%         fprintf('\n Downrange distance        = %4.2f [km]',x2(end))
%         fprintf('\n Drag loss                 = %4.2f [km/s]',vD2(end))
%         fprintf('\n Gravity loss              = %4.2f [km/s]',vG2(end))
%         fprintf('\n');
    end
        
    tf = t0 + tburn1;  % The time when propellant is completely burned
    %and the thrust goes to zero
    t_range = [t0,tf];  % Integration interval
     
%     figure(4)
%     plot(t,h2,'b');
%     hold on;
%     grid on;
%     xlabel('time[s]');
%     ylabel('Altitude[km]');
 
%     figure(5)
%     plot(t,x2,'b');
%     hold on;
%     grid on;
%     xlabel('time[s]');
%     ylabel('Downrange Distance[km]');
    
%     figure(6)
%     plot(t,gamma2,'b');
%     hold on;
%     grid on;
%     xlabel('time[s]');
%     ylabel('gamma [radians]');
end


