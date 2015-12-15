function dfdt = RocketDynEq(t,y)
global m0 g0 T A Cd rh0 H0 Re hgr_turn md
v  =  y(1);     % Velocity
gm =  y(2);     % Flight path angle
x  =  y(3);     % Downrange distance
h  =  y(4);     % Altitude
vD =  y(5);     % Velocity loss due to drag
vG =  y(6);     % Velocity loss due to gravity
% Equations of motion of a gravity turn trajectory
      m = m0 - md*t;  % Vehicle mass
% else
%     m = mf;          % Burnout mass
%     T = 0;           % No more thrust is generated
% end
g  = g0/(1 + h/Re)^2;          % Gravitational variation with altitude
rh = rh0*exp(-h/H0);            % Atmospheric density exponential model
D = 1/2 * rh*v^2 * A * Cd;      % Drag force

% Rocket starts the gravity turn when h = hgr_turn
if h <= hgr_turn % Vertical flight
    dv_dt  = T/m - D/m - g;           % (Q3)
    dgm_dt = 0;                       % (Q4)
    dx_dt  = 0;                       % (Q1)
    dh_dt  = v;                       % (Q2)
    dvG_dt = -g;                      % Gravity loss rate [m/s^2]
else                                  % Gravity turn
    dv_dt  = T/m - D/m - g*sin(gm); %(Q3)
    dgm_dt = -1/v*(g - v^2/(Re + h))*cos(gm); %(Q4)
    dx_dt  = Re/(Re + h)*v*cos(gm) + 463*sin(gm);  %(Q1)
    dh_dt  = v*sin(gm) + 463*cos(gm); % (Q2) Adding earth's rotation speed
    dvG_dt = -g*sin(gm);              % Gravity loss rate [m/s^2]
end
    dvD_dt = -D/m;           % Drag loss rate [m/s^2]
dfdt = [ dv_dt,dgm_dt, dx_dt,dh_dt, dvD_dt, dvG_dt]';
return