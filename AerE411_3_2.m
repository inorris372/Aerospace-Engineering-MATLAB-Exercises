%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Ian Norris           %
%     Propulsion Assignment 3    %
%           Problem 2            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 1.4;
Tt2 = 918.7;
Tt1 = 500;
taub = Tt2/Tt1;
upL = zeros(401);
M2 = zeros(401);
piB = zeros(401);


for M1 = 0:0.001:.4
    upL(int32(M1*1000+1))= M1^2*(1+(gamma-1)/2*M1^2)/(1+gamma*M1^2)^2*taub;
    M2(int32(M1*1000+1)) = sqrt((2*gamma*upL(int32(M1*1000+1))-1+...
        sqrt((2*gamma*upL(int32(M1*1000+1))-1)^2+4*upL(int32(M1*1000+1....
        ))*((gamma-1)/2-gamma^2*upL(int32(M1*1000+1)))))/(2*((gamma-1.....
        )/2-gamma^2*upL(int32(M1*1000+1)))));
    piB(int32(M1*1000+1)) = (1+gamma*M1^2)/(1+gamma*M2(int32(M1*1000+1)...
        )^2)*((1+(gamma-1)/2*M2(int32(M1*1000+1))^2)/(1+(gamma-1)/2*M1....
        ^2))^(gamma/(gamma-1));
end
M1 = (0:0.001:.4);

figure(1)
plot(M1,M2)
axis([0 0.4 0 1])
title('M1 vs M2')
xlabel('M1')
ylabel('M2')

figure(2)
plot(M1,piB)
axis([0 0.4 .8 1])
title('M1 vs piB')
xlabel('M1')
ylabel('piB = Pt2/Pt1')