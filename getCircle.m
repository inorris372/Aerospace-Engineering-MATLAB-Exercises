function [x,y] = getCircle(c1,c2, r)
numPoints=1000; %Number of points making up the circle
i = linspace(0,2*pi,numPoints);   %Angle Array
radius=r;      %Radius of the circle

%Define circle in polar coordinates (angle and radius)
rho=ones(1,numPoints)*radius; %Radius should be 1 for all 100 points

%Convert polar coordinates to Cartesian for plotting
[a,b] = pol2cart(i,rho); 
x = a - c1;
y = b + c2;

end
