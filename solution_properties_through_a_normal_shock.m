%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         solution properties through a normal shock                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


prompt = 'What is the value of M1?';
result = input(prompt);
z = result;
%Mach 1 speed value is set by user
gamma = 1.4;

M1 = 1:.01:(2*z);

M2 = sqrt((M1.^2+2/(gamma-1))./((2*gamma/(gamma-1)).*M1.^2-1));
%Mach speed at Position 2
P2P1rat = (1+gamma.*M1.^2)./(1+gamma.*M2.^2);
%Ratio of Pressure 1 to Pressure 2
T2T1rat = (1+(gamma-1)./2.*M1.^2)./(1+(gamma-1)./2.*M2.^2);
%Ratio of Temperature 1 to Temperature 2
p2p1rat = (M1./M2).*((1+((gamma-1)./2).*M2.^2)./(1+((gamma-1)./2).*M1.^2));
%Ratio of Density 1 to Density 2
Pt2Pt1rat = (1+gamma*M1.^2)./(1+gamma.*M2.^2).*((1+(gamma-1)./2.*M2.^2) ...
    ./(1+(gamma-1)./2.*M1.^2)).^(gamma./(gamma-1));

    
figure(1)
plot(M1,M2,'y')
hold on
plot(M1,P2P1rat,'m')
hold on
plot(M1,T2T1rat,'r')
hold on
plot(M1,p2p1rat,'c')
hold on
plot(M1,Pt2Pt1rat,'b')
axis([1 3 0 5])
legend('M2','P2/P1','T2/T1','p2/p1','Pt2/Pt1')
