% This calculate Pressure drop in pipeline...
% It adopt Colebroke Method to calculate friction factor in nonlaminar condition 
% using Newton Raphson Numerical Method of nonlinear algebraic equation.
% Reference book: Introduction to Numerical Method in Chemical Engineering

clc;clear;
%% Fluid Properties and condition
% Fluid=Air, 
T=25;           %degree Celcius
P=1;            %in atm
rho=1.23;       %fluid density in kg/m2 @T=25oC and 1atm
meu=1.79e-5;    %fluid viscosity in kg/m-s @T=25oC and 1atm

%% Pipe Properties
D=4e-3;            %in meter
Vdash=50;       %in m/s
epsilon=0.0015e-3;  %pipe roughness in m
% L=1;            %pipe length

% Calculating Raynold number
Re=(rho*Vdash*D)/meu;
if Re<=2100
    ff=64/Re;
else
    
    ff=0.001;       %initial guess for frictional factor
    
    fterm1=(epsilon/D)/3.7;
    fterm2=2.51/Re;
    fprimet1=2*fterm2;
    ffnew=0.01;
    
    iter=0;
    while abs(ffnew-ff)>1e-6
        ff=ffnew;
        f=(1/sqrt(ff))+2*log10(fterm1+(fterm2/sqrt(ff)));
        
        fprime=-0.5*ff^-1.5+((fprimet1*(-0.5*ff^-1.5))/(fterm1+fterm2*ff^-0.5));
        ffnew=ff-(f/fprime);
        
        iter=iter+1;
    end
end

%% Calculating Pressure Drop over a pipe length of 10m
L=1:1:10;
DeltaP=(ff*L*(Vdash^2)*rho)/(2*D);

plot(L,DeltaP)
xlabel('Pipe Length (meters)');
ylabel('Pressure Drop (N/m2)');
title('Pressure Drop against Pipe Length')
