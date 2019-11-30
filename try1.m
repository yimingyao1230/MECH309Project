%%%                 MECH 309 - Numerical Methods in Mech Eng                 %%%
 
% Presented to Prof Siva Nadarajah Winter 2019 - November 22th
 
%Yiming Yao 260769906
%Zechen Ren 260765431
 
clear all
close all
clc
 
%% Known Variables
 
gamma = 1.4; %specific heat ratio for air
R = 287.058; %J*kg^1*K^1 gas constant
Tinf = 293; %K freestream static temperature
Pinf = 100; %kN/m^2 freestream static pressure
C = 340; %m/s speed of sound, assumed constant ???
x = 50; %x-direction domain
y = 50; %y-direction domain
Minf = 0.2; %Mach number of freestream
Uinf = Minf*sqrt(gamma*R*Tinf); %Flow speed of freestream 
toc = 0.08; % thickness ratio
xle  = 20;
xte  = 21;
%% Initializing Variables 
dx = 0.1; % grid discrete distance
dy = 0.1; % grid discrete distance

Nx = x/dx; % x direction grid
Ny = y/dy; % y direction grid
 
phi = zeros (Nx*Ny,1); % phi initialization
phiNew = zeros (Nx*Ny,1); % phi initialization
miu = zeros (Nx*Ny,1); % miu initialization
cp = zeros (Nx*Ny,1); % cp initialization
plotphi = zeros(Nx,Ny);

a = zeros (Nx*Ny,1);
b = zeros (Nx*Ny,1);
c = zeros (Nx*Ny,1);
d = zeros (Nx*Ny,1);
e = zeros (Nx*Ny,1);
f = zeros (Nx*Ny,1);
g = zeros (Nx*Ny,1);
 
j = -1; % y direction tracing
i = -1; % x direction tracing
 
u_ = -1; % local phi derivative to x initialization
v_ = -1; % local phi derivative to y initialization
 
m = -1; % Mech number locally
A = zeros (Nx*Ny,1); % initialize A matrix
 
error = Inf;
tol = 1E-4;
count = 0;

xspan = linspace(0,x,Nx);
dydx = toc * (-4 * xspan + 82);
dydx(xspan<xle | xspan>xte ) = 0;

errorlist = nan(1,1000);

%% solve
while error > tol
    %lower boundary
    for i = 1:Nx
        phi(i,1) = phi(i+Nx,1) - Uinf * dydx(i) * dy;
    end

    % Calculate A & miu
    for j = 2:(Ny-1)
        for i = 2:(Nx-1)
            loc = (j-1) * Ny + i;
            A(loc,1) = (1 - Minf^2) - (gamma + 1) * Minf^2 / Uinf * (phi(loc+1,1) - phi(loc-1,1) ) / ( 2 * dx);
            if A(loc,1) > 0
                miu(loc,1) = 0;
            else
                miu(loc,1) = 1;
            end
        end
    end
    phiNew = phi;
    % factor
    for j = 2:(Ny-1)
        for i = 3:(Nx-1)
            loc = (j-1) * Ny + i;
            c(loc,1) = 1 / dy^2;
            g(loc,1) = miu(loc-1,1) * A(loc-1,1) / dx^2;
            d(loc,1) = (1-miu(loc,1)) * A(loc,1) / (dx)^2 - 2 * miu(loc-1,1) * A(loc-1,1) / (dx)^2;
            e(loc,1) = (1-miu(loc,1)) * A(loc,1) / (dx)^2;
            b(loc,1) = 1 / dy^2;
            a(loc,1) = miu(loc-1,1) * A(loc-1,1) / (dx)^2 - 2 * ((1-miu(loc,1))*A(loc,1))/(dx)^2 - 2/(dy)^2;
            phi(loc,1) = -(c(loc,1) * phi(loc-Nx,1) ...
                + g(loc,1) * phi(loc-2,1) ...
                + d(loc,1) * phi(loc-1,1) ...
                + e(loc,1) * phi(loc+1,1) ...
                + b(loc,1) * phi(loc+Nx,1) ...
                ) ./ a(loc,1); 
        end
    end

    error = max(max(abs(phiNew - phi)));
    count = count + 1;
    errorlist (count) = error;
end
%% plot
for j = 1 : Ny
    for i = 1 : Nx
        loc = (j-1) * Ny + i;
        plotphi(i,j) = phi(loc,1);
    end
end
pcolor(plotphi)
colorbar
% semilog (errorlist);


