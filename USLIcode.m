%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Euler Angle USLI Code            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USLI Rocket Team Euler Angle Code Information

%Created by:  Ian Norris
%Date created: 4/7/2014
%How to operate?   First, hard code in name of text file that will be used
%                  to read in the initial body frame accelerations and 
%                  angles.  The program then will read this data, and will
%                  hopefully convert all data into Euler Angle readouts
%                  and should output graphs of Euler Angle Readings vs time
%                  , altitude, etc...
filename = 'Uslitest.txt';
delimeterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimeterIn,headerlinesIn);
[b,c] = size(A.data);
%Body Frame Angle Matrices Derived from Text File are associated with 
%appropriate variables in the matrix below.
for i = 1:c/b
    for j = 1:b
        if j==1
            Ux(i) = A.data(j,i);
        elseif (j == 2)
            Uy(i) = A.data(j,i);
        elseif (j == 3)
            Uz(i) = A.data(j,i);
        elseif (j == 4)
            dP(i) = A.data(j,i);
        elseif (j == 5)
            dQ(i) = A.data(j,i);
        elseif (j == 6)
            dR(i) = A.data(j,i);
        else
            Time(i) = A.data(j,i);
        end
        if (i > 1)
            if (A.data(j,i)~=A.data(j,i-1))
                counter = 1;
            else
                counter = 0;
            end
        end
    end
    while(counter == 1)
        total = total+A.data(j,i);
    end
    
    P(i) = sum(P)+Time(i)*dP(i);  %psi body frame angle
    Q(i) = sum(Q)+Time(i)*dQ(i);  %theta body frame angle
    R(i) = sum(R)+Time(i)*dR(i);  %phi body frame angle
% The Rotation matrix values are calculated below for each iteration.
    
Rotation(i) = [cos(Q(i))*cos(R(i)) sin(P(i))*sin(Q(i))*cos(R(i))...
        -cos(P(i))*sin(R(i)) cos(P(i))*sin(Q(i))*cos(R(i))+sin(P(i))*...
        sin(R(i)); cos(Q(i))*sin(R(i)) sin(P(i))*sin(Q(i))*sin(R(i))+...
        cos(P(i))*cos(R(i)) cos(P(i))*sin(Q(i))*sin(R(i))-sin(P(i))*...
        cos(R(i)); -sin(Q(i)) sin(P(i))*cos(Q(i)) cos(P(i))*cos(Q(i))];
end




