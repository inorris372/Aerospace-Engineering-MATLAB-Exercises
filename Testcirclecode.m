numPoints=100; %Number of points making up the circle
radius=1;      %Radius of the circle

%Define circle in polar coordinates (angle and radius)
theta=linspace(0,2*pi,numPoints); %100 evenly spaced points between 0 and 2pi
rho=ones(1,numPoints)*radius; %Radius should be 1 for all 100 points

%Convert polar coordinates to Cartesian for plotting
[X,Y] = pol2cart(theta,rho); 

%Plot a red circle
plot(X,Y,'r-','linewidth',2);
axis square


%Altering the radius allows us to plot a second circle within the first one
radius=0.5;
rho=ones(1,numPoints)*radius;
[X,Y] = pol2cart(theta,rho);
hold on
plot(X,Y,'k-','linewidth',2);
hold off