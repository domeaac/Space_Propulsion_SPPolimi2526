clear
clc
close all

% Add necessary paths (user-defined function)
SetPath();

D_in=8;
D_out=8;
Dt=2;
alpha_conv=deg2rad(45);
alpha_div=deg2rad(25);
a=1;
b=1;
RaoL=2;

Length=6;
Pt=3e6;
Pc=9e6;

S=170e6;
S_conv=170e6;

rho=8000;
rho_conv = 8000; 

nozzleGeometry = [D_in,D_out,Dt,alpha_conv,alpha_div,a,b,RaoL];

[m_tot, D_in_esterno,t]=SizingStructuralRocket(Pc,Pt,nozzleGeometry,Length,S,rho,S_conv,rho_conv)