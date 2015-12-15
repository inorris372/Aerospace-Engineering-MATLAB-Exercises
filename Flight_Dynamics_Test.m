%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Flight Dynamics Graphs        %
%             Ian Norris             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C=1;
D=1000000;
%t=0;
alpha = exp(-.00242*t)*C*cos(1.5438*t+D);
plot(t,alpha)
