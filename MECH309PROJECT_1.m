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
c = 340; %m/s speed of sound, assumed constant ???
x = 50; %x-direction domain
y = 50; %y-direction domain
Minf = 0.2; %Mach number of freestream
Uinf = Minf * c; %Flow speed of freestream 

%Initializing Variables 
Nx = 60; % x direction grid
Ny = 60; % y direction grid
 
dx = x/Nx; % grid discrete distance
dy = y/Ny; % grid discrete distance
 
phi = zeros (Nx,Ny); % phi initialization
miu = zeros (Nx,Ny); % miu initialization
 
j = -1; % y direction tracing
i = -1; % x direction tracing
 
u_ = -1; % local phi derivative to x initialization
v_ = -1; % local phi derivative to y initialization

m = -1; % Mech number locally
A = zeros (Nx,Ny); % initialize A matrix

error = Inf;
count = 0;

% START
%Boundary conditon
while (error > 10^-2 && count < 1000)
    count = count + 1;
    phikm1 = phi; 

%update phi at second bottom row, neummann at x != 20,21
%update phi at second bottom row, neummann at x  = 20,21
for  k = 2:Nx-1 
    if k * x / Nx < 20 || k * x / Nx > 21
        phi ( Nx, k) = phi ( Nx -1 , k) ;
    else
        dydx = CalAirfoil(k * x / Nx);
        phi ( Nx, k) = phi ( Nx -1 , k) + dy *  Uinf * dydx;
    end
end

%A judgement

u_ = (phi(j,i-1) - phi(j,i+1))/(2*dx) ;
v_ = (phi(j-1,i) - phi(j-1,i))/(2*dy) ;
u = (u_ + Uinf);
U = sqrt (u^2 + v_^2);
m = U / c;
if m > 1 % supersonic locally
    A (j,i) = (1 - Minf)^ 2 - (gamma + 1) * Minf^2 / Uinf * (phi (j,i) - phi (j,i-2) ) / ( 2* dx);
else % subsonic flow
    A (j,i) = (1 - Minf)^ 2 - (gamma + 1) * Minf^2 / Uinf * (phi (j,i+1) - phi (j,i-1) ) / ( 2* dx);
end
 
if A(j,i) > 0
    miu(j,i) = 0;
else
    miu(j,i)= 1;
end
    

%Point Gaussi-Seidal
for i = 3 : Nx-1
    for j = 2 : Ny-1
        phi (j,i) = -(...
             ( ( 1-miu(j,i) * A(j,i))/(dx)^2 ) * phi(j,i+1)...
            +( ( miu(j,i-1) * A(j,i-1) )/(dx)^2 ) * phi(j,i-2) ...
            +( ( 1/(dy)^2) * phi(j+1,i) )...
            +( ( 1/(dy)^2) * phi(j-1,i) )...
            +( ( (1-miu(j,i)) * A(j,i) / (dx)^2 - 2*miu(j,i-1) * A(j,i-1) / (dx)^2) * phi(j,i-1))...
            ) / ( miu(j,i-1)*A(j,i-1) / (dx)^2 - 2*(1-miu(j,i))*A(j,i) / (dx)^2 - 2/(dy)^2 );
    end
end

error = norm (phikm1 - phi);
end
%Plots