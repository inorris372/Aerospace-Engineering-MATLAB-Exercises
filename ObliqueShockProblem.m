%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Oblique Shock Problem            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part A
%Input Parameters should be the shock angle and the upstream flow
%conditions.  Output parameters should be the downstream Mach number and
%the deflection angle.  This will output a plot.
clear
clc
usigma = input('What is the value of the shock angle to plot/guess? \n');
%1st Mach speed is set in this problem
M1a = 2;
M1b = 1.5;
gamma = 1.4;
sigma = asin(1./M1a):.01:pi/2;
sigma2 = asin(1./M1b):.01:pi/2;

turn1 = atan(2./tan(sigma).*(((M1a.^2).*sin(sigma).^2-1)./...
    ((M1a.^2).*(gamma+cos(2.*sigma))+2)));
turn2 = atan(2./tan(sigma2).*(((M1b.^2).*sin(sigma2).^2-1)./...
    ((M1b.^2).*(gamma+cos(2.*sigma2))+2)));
turn1A = 2.*cot(sigma).*(M1a.^2.*sin(sigma).^2-10./(M1a.^2.*(gamma+...
    cos(2.*sigma))+2))-tan(turn1);
turn2A = 2.*cot(sigma2).*(M1b.^2.*sin(sigma2).^2-10./(M1b.^2.*(gamma+...
    cos(2.*sigma2))+2))-tan(turn2);

S1 = radtodeg(sigma);
S2 = radtodeg(sigma2);
T1 = radtodeg(turn1);
T2 = radtodeg(turn2);


figure(1)
plot(T1,S1,'r')
title('Turning Angle vs. Shock Angle')
xlabel('Delta (deg)')
ylabel('Sigma (deg)')
axis([0 30 20 100])
hold on
plot(T2,S2,'b')                         
legend('1st M1 value','2nd M1 value')

%Write a program to input a specified turning angle, upstream Mach number
%and flow conditions and compute the shock angle and downstream flow
%conditions. Your code should take the solution branch (i.e. weak or
%strong) as an input parameter.

upMach = input('What is the value of the upstream Mach Number? \n');
d1 = input('What is the the turning angle guess @ M 1.5? \n');
Shock1 = input('Is this the Weak or Strong Branch? (1:Weak, 2:Strong)\n');
d2 = input('What is the the turning angle guess @ M 2? \n');
Shock2 = input('Is this the Weak or Strong Branch? (1:Weak, 2:Strong)\n');

turn1prime = diff(sigma)./diff(turn1A);
turn2prime = diff(sigma2)./diff(turn2A);
% midsigmaA = radtodeg(max(turn1));
% midsigmaB = radtodeg(max(turn2));



d1A = turn1A(1:end-1);
d2A = turn2A(1:end-1);

for x = 1:2
    if x == 1 
        for z = 1:size(d1A)
           [sigmaPA, Error1] = Newton_Function(d1A(z),turn1prime(z)...
               ,usigma,Shock1,.0001);
        end
    else
        for y = 1:size(d2A)
           [sigmaPB, Error2] = Newton_Function(d2A(y),turn2prime(y)...
               ,usigma,Shock2,.0001);
        end
    end
end

plot(d1,sigmaPA,'o','MarkerSize',10);
hold on
plot(d2,sigmaPB,'o','MarkerSize',10);
     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        




