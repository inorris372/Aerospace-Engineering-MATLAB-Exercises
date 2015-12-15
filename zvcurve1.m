% zvcurve1.m      May 5, 2008

% plots zero relative velocity curves
% through the equilibrium points of the
% circular-restricted Earth-Moon
% three body problem

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

clc; home;

fprintf('\n             program zvcurve1\n');

fprintf('\ngraphics display of zero relative velocity');
fprintf('\ncurves through the crtbp equilibrium points\n\n');

while(1)
    fprintf('\n <1> L1 contour\n\n');

    fprintf(' <2> L2 contour\n\n');

    fprintf(' <3> L3 contour\n\n');

    fprintf(' <4> L4 and L5 contour\n\n');

    fprintf(' <5> all points contour\n\n');

    fprintf(' selection (1, 2, 3, 4 or 5)\n\n');

    icflg = input('? ');

    if (icflg >= 1 && icflg <= 5)
        break;
    end
end

fprintf('\n\n  working ...\n');

% define Earth-Moon mass ratio

mu = 0.012155099;

% create x and y mesh data points

xmin = -1.5;

xmax = +1.5;

ymin = -1.5;

ymax = +1.5;

delta = 0.005;

[x, y] = meshgrid(xmin: delta: xmax, ymin: delta: ymax);

% create z data

x1 = - mu;

x2 = 1 - mu;

for i = 1: 1: size(x)
    for j = 1: 1: size(x)

        r1sqr = (x(i, j) - x1)^2 + y(i, j)^2;

        r2sqr = (x(i, j) - x2)^2 + y(i, j)^2;

        z(i, j) = (1 - mu) * (r1sqr + 2 / sqrt(r1sqr)) ...
            + mu * (r2sqr + 2 / sqrt(r2sqr)) - mu * (1 - mu);
    end
end

% define contour levels

% L1 libration point

v(1) = 2 * 1.5941913671;

% L2 libration point

v(2) = 2 * 1.5860980401;

% L3 libration point

v(3) = 2 * 1.5060758307;

% L4 and L5 libration points

v(4) = 2 * 1.4939963237;

v(5) = 2 * 1.4939963237;

v(4) = 3.0;

% plot zero relative-velocity contours

% clevel = [v(3) v(3)];

switch icflg
    case 1
        clevel = [v(1) v(1)];
    case 2
        clevel = [v(2) v(2)];
    case 3
        clevel = [v(3) v(3)];

        [c, h] = contour(x, y, z, v(3));

        clabel(c, v(3));
    case 4
        clevel = [v(4) v(5)];

        [c, h] = contour(x, y, z, v(4));

        clabel(c, v);

    case 5
        clevel = [v(1) v(2) v(3) v(4)];

        [c, h] = contour(x, y, z, v);

        clabel(c, v);
end

%contour (x, y, z, clevel);

axis square;

axis ([xmin xmax ymin ymax]);

grid off;

hold on;

% plot body locations

% Earth

plot(-mu, 0, '.b');

% Moon

plot(1 - mu, 0, '.b');

% plot location of libration points

switch icflg
    case 1
        % L1

        xl = 0.83689291982029;

        yl = 0;

        plot(xl, yl, '*r');

        title('Zero Relative Velocity Curve Through L1', 'FontSize', 16);
    case 2
        % L2

        xl = 1.15569952178709;

        yl = 0;

        plot(xl, yl, '*r');

        title('Zero Relative Velocity Curve Through L2', 'FontSize', 16);
    case 3
        % L3

        xl = -1.00506452627984;

        yl = 0;

        plot(xl, yl, '*r');

        title('Zero Relative Velocity Curve Through L3', 'FontSize', 16);
    case 4
        % L4 and L5

        xl = 0.4878449010;

        yl = 0.86602540378444;

        plot(xl, yl, '*r');

        xl = 0.4878449010;

        yl = -0.86602540378444;

        plot(xl, yl, '*r');

        title('Zero Relative Velocity Curve through L4 and L5', 'FontSize', 16);
    case 5
        xl = 0.83689291982029;

        yl = 0;

        plot(xl, yl, '*r');

        xl = 1.15569952178709;

        yl = 0;

        plot(xl, yl, '*r');

        xl = -1.00506452627984;

        yl = 0;

        plot(xl, yl, '*r');

        xl = 0.4878449010;

        yl = 0.86602540378444;

        plot(xl, yl, '*r');

        xl = 0.4878449010;

        yl = -0.86602540378444;

        plot(xl, yl, '*r');

        title('Zero Relative Velocity Curves through L1-L5', 'FontSize', 16);
end

xlabel('x coordinate');

ylabel('y coordinate');

print -depsc -tiff -r300 zvcurve1.eps

