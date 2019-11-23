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
toc = 0.08; % thickness ratio
syms x_; 
yx = toc * (-2 * x_^2 + 82*x_ -840); % circular-arc defination for airfoil
dydx = diff(yx,x);


%Initializing Variables 
Nx = 60;
Ny = 60;

dx = x/Nx;
dy = y/Ny;

phi = zeros (Nx,Ny);
miu = zeros (Nx,Ny);

j = -1;
i = -1;

error = Inf;
count = 0;
%Boundary conditon
while (error > 10^-2)
    count = count + 1;
%phi at third column (i = 3) BC enforced at intialization 
%phi at second top row (j = Ny-1) BC enforced at intialization 
%phi at second right column (i = Nx-1) BC enforced at intialization 

%update phi at second bottom row, neummann at x != 20,21
%update phi at second bottom row, neummann at x  = 20,21
for  k = 2:Nx-1 
    if k * x / Nx < 20 || k * x / Nx > 21
        phi ( Nx, k) = phi ( Nx -1 , k) ;
    else
        phi ( Nx, k) = phi ( Nx -1 , k) + dy *  Uinf * subs(dydx,x_,(k * x / Nx));
    end
end
%A judgement

u = diff (phi(j,i),i);

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

error = norm (phi);
end
%Plots