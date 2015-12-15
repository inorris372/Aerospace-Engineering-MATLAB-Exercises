%Zach Willnerd
%AerE 451
%Homework 4
%Problem 2 (Chobotov Prob 5.2)

clear all
clc

r=6378; %km
rA=9*r; %km
rB=16*r; %km
rC=16*r; %km
rD=25*r; %km
rE=25*r; %km
mu=3.986*10^5; %km^3/s^2

%% Double Hohmann Transfer (A-B-C-D)

Vc_A=sqrt(mu/rA); %km/s
Vc_B=sqrt(mu/rB); %km/s
Vc_C=sqrt(mu/rC); %km/s
Vc_D=sqrt(mu/rD); %km/s
Vc_E=sqrt(mu/rE); %km/s

Vt_AB1=sqrt((2*mu/rA)-((2*mu)/(rA+rB))); %km/s
Vt_AB2=sqrt((2*mu/rB)-((2*mu)/(rA+rB))); %km/s
Vt_CD1=sqrt((2*mu/rC)-((2*mu)/(rC+rD))); %km/s
Vt_CD2=sqrt((2*mu/rD)-((2*mu)/(rC+rD))); %km/s

DeltaVA=Vt_AB1-Vc_A; %km/s
DeltaVB=Vc_B-Vt_AB2; %km/s
DeltaVC=Vt_CD1-Vc_C; %km/s
DeltaVD=Vc_D-Vt_CD2; %km/s

DeltaV_totalkm=DeltaVA+DeltaVB+DeltaVC+DeltaVD; %km/s
DeltaV_totalm=DeltaV_totalkm*1000; %m/s

fprintf('The total velocity for a double Hohmann Transfer (A-B-C-D) = %.15f m/s \n\n',DeltaV_totalm);

%% Single Hohmann Transfer (A-D)

Vc_A=sqrt(mu/rA); %km/s
Vc_D=sqrt(mu/rD); %km/s

Vt_AD1=sqrt((2*mu/rA)-((2*mu)/(rA+rD))); %km/s
Vt_AD2=sqrt((2*mu/rD)-((2*mu)/(rA+rD))); %km/s

DeltaVA=Vt_AD1-Vc_A; %km/s
DeltaVD=Vc_D-Vt_AD2; %km/s

DeltaV_totalkm=DeltaVA+DeltaVD; %km/s
DeltaV_totalm=DeltaV_totalkm*1000; %m/s

fprintf('The total velocity for a Single Homhman Transfer (A-D) = %.15f m/s \n\n',DeltaV_totalm);

%% Bi-Elliptic Transfer (A-B-E)

Vc_A=sqrt(mu/rA); %km/s
Vc_E=sqrt(mu/rE); %km/s

Vt_AB1=sqrt((2*mu/rA)-((2*mu)/(rA+rB))); %km/s
Vt_AB2=sqrt((2*mu/rB)-((2*mu)/(rA+rB))); %km/s
Vt_BE1=sqrt((2*mu/rB)-((2*mu)/(rB+rE))); %km/s
Vt_BE2=sqrt((2*mu/rE)-((2*mu)/(rB+rE))); %km/s

DeltaVA=Vt_AB1-Vc_A; %km/s
DeltaVB=Vt_BE1-Vt_AB2; %km/s
DeltaVE=Vc_E-Vt_BE2; %km/s

DeltaV_totalkm=DeltaVA+DeltaVB+DeltaVE; %km/s
DeltaV_totalm=DeltaV_totalkm*1000; %m/s

fprintf('The total velocity for a Bi-Elliptic Transfer (A-B-E) = %.15f m/s \n\n',DeltaV_totalm);
