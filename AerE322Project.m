%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AerE 322 Project - Ian Norris    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Member Stiffnes in global Coordinates
A = 6; %Cross Sectional Area
I = 18;  %Moment of Inertia
E = 29000; %Young's Modulus for both members

for i = 1:2
    if i==1
        L = 300;
        theta = acos(.6);
        sub1 = A*L^2/I*cos(theta);
        sub2 = 12*sin(theta);
        sub3 = sub1*tan(theta);
        sub4 = sub2*cot(theta);
        K1 = E.*I./(L.^3).*[(sub1*cos(theta)+sub2*sin(theta)),...
    (sub1*sin(theta)-sub2*cos(theta)),(-L/2*sub2),...
    -(sub1*cos(theta)+sub2*sin(theta)),...
    (-sub1*sin(theta)+sub2*cos(theta)),(-L/2*sub2); 
    (sub1*sin(theta)-sub2*cos(theta)),(sub3*sin(theta)+sub4*cos(theta)),...
    (L/2*sub4),(-sub1*sin(theta)+sub2*cos(theta)),...
    -(sub3*sin(theta)+sub4*cos(theta)),(L/2*sub4); (-L/2*sub2),...
    (L/2*sub4),(4*L^2),(L/2*sub2),(-L/2*sub4),(2*L^2); 
    (-sub1*cos(theta)-sub2*sin(theta)),...
    (-sub1*sin(theta)+sub2*cos(theta)),(L/2*sub2),...
    (sub1*cos(theta)+sub2*sin(theta)),(sub1*sin(theta)-sub2*cos(theta)),...
    (L/2*sub2); (-sub1*sin(theta)+sub2*cos(theta)),...
    -(sub3*sin(theta)+sub4*cos(theta)),(-L/2*sub4),...
    (sub1*sin(theta)-sub2*cos(theta)),(sub3*sin(theta)+sub4*cos(theta)),...
    (-L/2*sub4); (-L/2*sub2),(L/2*sub4),(2*L^2),(L/2*sub2),(-L/2*sub4),...
    4*L^2];

    M1 = E*I/(L^3)*[A*L^2/I 0 0 -A*L^2/I 0 0; 0 12 6*L 0 -12 6*L;
        0 6*L 4*L^2 0 -6*L 2*L^2; -A*L^2/I 0 0 A*L^2/I 0 0;
        0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2];
    
    N1 = ((E*I)/L^3)*M1
    
    T1 = [cos(theta) sin(theta) 0 0 0 0; -sin(theta) cos(theta) 0 0 0 0;
    0 0 1 0 0 0; 0 0 0 cos(theta) sin(theta) 0;
    0 0 0 -sin(theta) cos(theta) 0; 0 0 0 0 0 1]
    else
        L = 216.33;
        theta = -acos(.832);
        sub1 = A*L^2/I*cos(theta);
        sub2 = 12*sin(theta);
        sub3 = sub1*tan(theta);
        sub4 = sub2*cot(theta);
        K2 = E.*I./(L.^3).*[(sub1*cos(theta)+sub2*sin(theta)),...
    (sub1*sin(theta)-sub2*cos(theta)),(-L/2*sub2),...
    -(sub1*cos(theta)+sub2*sin(theta)),...
    (-sub1*sin(theta)+sub2*cos(theta)),(-L/2*sub2); 
    (sub1*sin(theta)-sub2*cos(theta)),(sub3*sin(theta)+sub4*cos(theta)),...
    (L/2*sub4),(-sub1*sin(theta)+sub2*cos(theta)),...
    -(sub3*sin(theta)+sub4*cos(theta)),(L/2*sub4); (-L/2*sub2),...
    (L/2*sub4),(4*L^2),(L/2*sub2),(-L/2*sub4),(2*L^2); 
    (-sub1*cos(theta)-sub2*sin(theta)),...
    (-sub1*sin(theta)+sub2*cos(theta)),(L/2*sub2),...
    (sub1*cos(theta)+sub2*sin(theta)),(sub1*sin(theta)-sub2*cos(theta)),...
    (L/2*sub2); (-sub1*sin(theta)+sub2*cos(theta)),...
    -(sub3*sin(theta)+sub4*cos(theta)),(-L/2*sub4),...
    (sub1*sin(theta)-sub2*cos(theta)),(sub3*sin(theta)+sub4*cos(theta)),...
    (-L/2*sub4); (-L/2*sub2),(L/2*sub4),(2*L^2),(L/2*sub2),(-L/2*sub4),...
    4*L^2];

    M2 = E*I/(L^3)*[A*L^2/I 0 0 -A*L^2/I 0 0; 0 12 6*L 0 -12 6*L;
        0 6*L 4*L^2 0 -6*L 2*L^2; -A*L^2/I 0 0 A*L^2/I 0 0;
        0 -12 -6*L 0 12 -6*L; 0 6*L 2*L^2 0 -6*L 4*L^2];
    
    N2 = ((E*I)/L^3)*M2
    
    T2 = [cos(theta) sin(theta) 0 0 0 0; -sin(theta) cos(theta) 0 0 0 0;
    0 0 1 0 0 0; 0 0 0 cos(theta) sin(theta) 0;
    0 0 0 -sin(theta) cos(theta) 0; 0 0 0 0 0 1]
    end
end
K1
K2

S = [K1(4,4)+K2(1,1),K1(4,5)+K2(1,2),K1(4,6)+K2(1,3); K1(5,4)+K2(2,1),... 
    K1(5,5)+K2(2,2),K1(5,6)+K2(2,3); K1(6,4)+K2(3,1),K1(6,5)+K2(3,2),...
    K1(6,6)+K2(3,3)]

Pf = [0 0 0]
P = [-150 100 120]
d = (P-Pf)/S

reshape(d,3,1)
%% Frame Disassembly in Global Coordinates

v1=[0;0;0;d(1,1);d(1,2);d(1,3)];
v2=[d(1,1);d(1,2);d(1,3);0;0;0];
u1 = T1*v1
u2 = T2*v2

Q1 = N1*u1
Q2 = N2*u2

sigma_a1 = Q1(1)/A
sigma_a2 = Q2(1)/A
epsilon_a1 = sigma_a1/E
epsilon_a2 = sigma_a2/E

F1 = transpose(T1)*Q1
F2 = transpose(T2)*Q2

c = 3;

X1 = -sqrt(20^2+10^2)*12:1:0;
X2 = 0:1:12*20;
sigmab1 = ((100*X1-120)*c)/I;
sigmab2 = ((100*X2-120)*c)/I;

figure(1)
plot(sigmab1)
xlabel('Length (inches)')
ylabel('Stress (ksi)')

figure(2)
plot(sigmab2)
xlabel('Length (inches)')
ylabel('Stress (ksi)')








