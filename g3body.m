% g3body.m      November 9, 2012

% graphics display of orbital motion in the
% circular-restricted Earth-Moon three body problem

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global mu

clc; home;

fprintf('\n             program g3body\n');

fprintf('\n< graphics display of three-body motion >\n\n');

while(1)
    
    fprintf('\n <1> periodic orbit about L1\n\n');

    fprintf(' <2> periodic orbit about L2\n\n');

    fprintf(' <3> periodic orbit about L3\n\n');

    fprintf(' <4> user input of initial conditions\n\n');

    fprintf(' selection (1, 2, 3 or 4)\n\n');

    icflg = input('? ');

    if (icflg >= 1 && icflg <= 4)
        
        break;
        
    end
    
end

switch icflg
    
    case 1
        
        % periodic l1 orbit (tr 32-1168, pp 25,29; 74)

        y(1) = 0.300000161;
        y(3) = 0;
        y(2) = 0;
        y(4) = -2.536145497;

        ti = 0;
        tf = 5.349501906;

        mu = 0.012155092;

        % set plot boundaries

        xmin = -1.5;
        xmax = +1.5;
        ymin = -1.5;
        ymax = +1.5;

    case 2
        
        % periodic l2 orbit (tr 32-1168, pp 31,34; 126)

        y(1) = 2.840829343;
        y(3) = 0;
        y(2) = 0;
        y(4) = -2.747640074;

        ti = 0;
        tf = 2 * 5.966659294;

        mu = 0.012155085;

        % set plot boundaries

        xmin = -3;
        xmax = +3;
        ymin = -3;
        ymax = +3;

    case 3
        
        % periodic l3 orbit (tr 32-1168, pp 37,39; 63)

        y(1) = -1.600000312;
        y(3) = 0;
        y(2) = 0;
        y(4) = 2.066174572;

        ti = 0;
        tf = 2 * 3.151928156;

        mu = 0.012155092;

        % set plot boundaries

        xmin = -2;
        xmax = +2;
        ymin = -2;
        ymax = +2;

    case 4
        
        % user input of initial conditions

        fprintf('\nplease input the x component of the radius vector\n');

        y(1) = input('? ');

        fprintf('\nplease input the y component of the radius vector\n');

        y(3) = input('? ');

        fprintf('\nplease input the x component of the velocity vector\n');

        y(2) = input('? ');

        fprintf('\nplease input the y component of the velocity vector\n');

        y(4) = input('? ');

        while(1)

            fprintf('\nplease input the value for the earth-moon mass ratio\n');

            mu = input('? ');

            if (mu > 0)
                break;
            end
        end

        while(1)

            fprintf('\nplease input the final time\n');

            tf = input('? ');

            if (abs(tf) > 0)
                break;
            end
        end

        % default values for plot boundaries

        xmin = -2;
        xmax = +2;
        ymin = -2;
        ymax = +2;

end

% request the integration step size

while(1)
    
    fprintf('\n\nplease input the integration step size\n');

    fprintf('(a value of 0.01 is recommended)\n');

    dt = input('? ');

    if (abs(dt) > 0)
        
        break;
        
    end
    
end

% set ode45 options

options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10);

% initialize

t2 = -dt;

npt = 0;

fprintf('\n\n  working ...\n');

while(1)

    t1 = t2;

    t2 = t1 + dt;

    [twrk, ysol] = ode45(@crtbp_eqm, [t1, t2], y, options);
    
    npt = npt + 1;

    xplot(npt) = ysol(length(twrk), 1);

    yplot(npt) = ysol(length(twrk), 3);

    y = ysol(length(twrk), 1:4);

    % check for end of simulation

    if (t2 >= tf)
        
        break;
        
    end
    
end

% plot trajectory

plot(xplot, yplot);

axis ([xmin xmax ymin ymax]);

axis square;

ylabel('y coordinate');

xlabel('x coordinate');

% label locations of Earth and Moon

hold on;

plot(-mu, 0, '*g');

plot(1 - mu, 0, '*b');

% label libration points

switch icflg
    
    case 1
        
        plot(0.836892919, 0, '.r');

        title('Periodic Orbit about the L1 Libration Point', 'FontSize', 16);
        
    case 2
        
        plot(1.115699521, 0, '.r');

        title('Periodic Orbit about the L2 Libration Point', 'FontSize', 16);
        
    case 3
        
        plot(-1.005064527, 0, '.r');

        title('Periodic Orbit about the L3 Libration Point', 'FontSize', 16);
        
    case 4
        
        title('User Defined Initial Conditions', 'FontSize', 16);
end

% create eps graphics file with tiff preview

print -depsc -tiff -r300 g3body.eps


