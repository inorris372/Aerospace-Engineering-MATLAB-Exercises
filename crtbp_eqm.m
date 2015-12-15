function ydot = crtbp_eqm (tsim, y)

% first order circular restricted
% three-body equations of motion

% required by g3body.m

% input

%  y(1) = x-component of position
%  y(2) = x-component of velocity
%  y(3) = y-component of position
%  y(4) = y-component of velocity

% output

%  ydot(1) = x component of velocity
%  ydot(2) = x component of acceleration
%  ydot(3) = y component of velocity
%  ydot(4) = y component of acceleration

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu

r1 = sqrt((y(1) + mu) ^ 2 + y(3) ^ 2);

r2 = sqrt((y(1) - 1 + mu) ^ 2 + y(3) ^ 2);

% integration vector

ydot = [y(2)
        2 * y(4) + y(1) - (1 - mu) * (y(1) + mu) / r1 ^ 3 - mu ...
    * (y(1) - 1 + mu) / r2 ^ 3
        y(4)
        -2 * y(2) + y(3) - (1 - mu) * y(3) / r1 ^ 3 - mu ...
    * y(3) / r2 ^ 3];
