%%%%%%%%%%%%%%%%%%%%%%%
%    AerE 331 Code    %
%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 2
%Part B

nume1 = 5;
denom1 = [1 1 25];
[A,B,C,D] = tf2ss(nume1,denom1);
A
B
C
D

%Part C
b1 = [5; 0];
b2 = [1; 0];
c1 = [0 1];
c2 = [0 5];
[nume2,denom2] = ss2tf(A,b1,c1,D);
trans1 = tf(nume2,denom2)
[nume3,denom3] = ss2tf(A,b2,c2,D);
trans2 = tf(nume3,denom3)

%Part D

%% Problem 3
%Part B
p2 = [-10+10i -10-10i];
Kp1 = place(A,b1,p2)
bsw = [1;0];
Kp2 = place(A,bsw,p2)

%Part C
trans3 = tf(nume2,denom2*8)
trans5 = tf(nume3,denom3*8)

%Part D
[a7,b7,c7,d7] = tf2ss(nume2,denom2);
sys2 = ss(a7,b7,c7,d7);
num4 = 4;
denom4 = [1 4 4];
trans4 = tf(num4,denom4);
[a3,b3,c3,d3] = tf2ss(num4,denom4);
sys1 = ss(a3,b3,c3,d3);
figure(7)
step(trans1)
hold on
u = feedback(sys2,sys1);
step(u)


%% Problem 4
%Part E

num4 = 4;
denom4 = [1 4 4];
trans4 = tf(num4,denom4);
figure(1)
step(trans4)
grid

x0 = [0 1];
[a3,b3,c3,d3] = tf2ss(num4,denom4);
sys1 = ss(a3,b3,c3,d3);
figure(2)
initial(sys1,x0)
grid
