%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Prandtl-Meyer Expansion      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computes exit Mach Number and pressure as a function of upstream Mach
%number and pressure (for a general set of input parameters).

theta1 = asin(1/M1);
mu1 = theta1;
C = -theta1 - sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1))...
    *sqrt(M1^2-1));
d = C + asin(1/M2)+sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)...
    /(gamma+1))*sqrt(M2^2-1));
Tt = T*(1+(gamma-1)/2*M^2);
Pt = P1*((1+gamma-1)/2*M1^2)^(gamma/(gamma-1));
P2 = Pt/(1+(gamma-1)/2*M2^2)^(gamma/(gamma-1));
