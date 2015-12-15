%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Ian Norris            %
%     Astro 451 Hw 5        %
%     Problem 5.2           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=6378; %km
mu = 3.986*10^5; %km^3/s^2
ra = 9*r; %km
rb = 16*r; %km
rc = 16*r; %km
rd = 25*r; %km
re = 25*r; %km

% Double Hohmann Transfer Section (A-B-C-D)

Vca = sqrt(mu/ra); %km/s
Vcb = sqrt(mu/rb); %km/s
Vcc = sqrt(mu/rc); %km/s
Vcd = sqrt(mu/rd); %km/s
Vce = sqrt(mu/re); %km/s
Vtab1 = sqrt((2*mu/ra)-((2*mu)/(ra+rb))); %km/s
Vtab2 = sqrt((2*mu/rb)-((2*mu)/(ra+rb))); %km/s
Vtcd1 = sqrt((2*mu/rc)-((2*mu)/(rc+rd))); %km/s
Vtcd2 = sqrt((2*mu/rd)-((2*mu)/(rc+rd))); %km/s
Deltava = Vtab1-Vca; %km/s
Deltavb = Vcb-Vtab2; %km/s
Deltavc = Vtcd1-Vcc; %km/s
Deltavd = Vcd-Vtcd2; %km/s
DeltaVtkm = Deltava+Deltavb+Deltavc+Deltavd; %km/s
DeltaVtm = DeltaVtkm*1000; %m/s
fprintf('Dbl. Hohmann Transfer (A-B-C-D) Vt = %.15f m/s\n\n',DeltaVtm);

% Single Hohmann Transfer (A-D)

Vca = sqrt(mu/ra); %km/s
Vcd = sqrt(mu/rd); %km/s
Vtad1 = sqrt((2*mu/ra)-((2*mu)/(ra+rd))); %km/s
Vtad2 = sqrt((2*mu/rd)-((2*mu)/(ra+rd))); %km/s
Deltava = Vtad1-Vca; %km/s
Deltavd = Vcd-Vtad2; %km/s
DeltaVtkm = Deltava+Deltavd; %km/s
DeltaVtm = DeltaVtkm*1000; %m/s
fprintf('Single Homhman Transfer (A-D) Vt = %.15f m/s \n\n',DeltaVtm);

% Bi-Elliptic Transfer (A-B-E)

Vca = sqrt(mu/ra); %km/s
Vce = sqrt(mu/re); %km/s
Vtab1 = sqrt((2*mu/ra)-((2*mu)/(ra+rb))); %km/s
Vtab2 = sqrt((2*mu/rb)-((2*mu)/(ra+rb))); %km/s
Vtbe1 = sqrt((2*mu/rb)-((2*mu)/(rb+re))); %km/s
Vtbe2 = sqrt((2*mu/re)-((2*mu)/(rb+re))); %km/s
Deltava = Vtab1-Vca; %km/s
Deltavb = Vtbe1-Vtab2; %km/s
Deltave = Vce-Vtbe2; %km/s
DeltaVtkm = Deltava+Deltavb+Deltave; %km/s
DeltaVtm = DeltaVtkm*1000; %m/s
fprintf('Bi-Elliptic Transfer (A-B-E) Vt = %.15f m/s \n\n',DeltaVtm);
