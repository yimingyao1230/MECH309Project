%%%                 MECH 309 - Numerical Methods in Mech Eng                 %%%

% Presented to Prof Siva Nadarajah Winter 2019 - November 22th

%Yiming Yao 260769906
%Zechen Ren 260765431
%

clear all
close all
clc
 
%Known Variables

gamma = 1.4; %specific heat ratio for air
R = 287.058; %J*kg^?1*K^?1 gas constant
Tinf = 293; %K freestream static temperature
Pinf = 100; %kN/m^2 freestream static pressure
x = 50; %x-direction domain
y = 50; %y-direction domain
Minf = 0.2; %Mach number of freestream
%Initializing Variables 
Nx = 60;
Ny = 60;

dx = x/NX;
dy = y/Ny;

phi = zeros (Nx,Ny);

%Boundary conditon

%Point Jacobi

%Plots