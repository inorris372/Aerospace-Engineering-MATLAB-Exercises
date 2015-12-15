%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Astrodynamics Project Code %
%         Ian Norris          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E1 = 0.11;  %Structure Ratio for Stage 1
E2 = 0.13;  %Structure Ratio for Stage 2

ISP1 = 280;  %ISP for Stage 1 Rocket
ISP2 = 455;  %ISP for Stage 2 Rocket

PiK = [(m2+3200)/m1,3200/m2];
x = [E1,E2];  %An array to store all Structure Ratio values
z = [ISP1,ISP2];  %An array to store all ISP values

-z(1).*.00981*ln(x(1)+(1-x(1))*PiK(1))-z(2).*.00981*ln(x(2)+(1-x(2))*PiK(2))=11.51752522;
% This equation is supposed to graph the masses of stage 1
% and stage 2, in order to optimize the rocket for best performance
% and lightest structural system.
plot(m1,m2)