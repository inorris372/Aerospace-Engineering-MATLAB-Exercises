function [ R ] = rstTOijk(x)
%Rotating S Frame to Earth Centered Inertial Frame 
%Given a position vector computes the velociety transformation
%See celestial mechanisc note set one.

normX=sqrt((x(1,1))^2+(x(2,1))^2+(x(3,1))^2);
partX=sqrt((x(1,1))^2+(x(2,1))^2);

R=[x(1,1)/normX,  (-1)*x(2,1)/partX, x(1,1)*x(2,1);
    x(2,1)/normX, x(1,1)/partX, x(2,1)*x(3,1);
    x(3,1)/normX, 0, (x(1,1)*x(1,1)+x(2,1)*x(2,1))];
end

