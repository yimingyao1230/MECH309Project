%%%                 MECH 309 - Numerical Methods in Mech Eng                 %%%
 
% Presented to Prof Siva Nadarajah Winter 2019 - November 22th
 
%Yiming Yao 260769906
%Zechen Ren 260765431
%
clear all
close all
clc
 
%% Known Variables
 
gamma = 1.4; %specific heat ratio for air
R = 287.058; %J*kg^1*K^1 gas constant
Tinf = 293; %K freestream static temperature
Pinf = 100000; %N/m^2 freestream static pressure
C = 340; %m/s speed of sound, assumed constant ???
x = 50; %x-direction domain
y = 50; %y-direction domain
Minf = 0.8; %Mach number of freestream
Uinf = Minf*sqrt(gamma*R*Tinf); %Flow speed of freestream 
toc = 0.08; % thickness ratio
xle  = 20; % Airfoil Leading Edge
xte  = 21; % Airfoil Trailing Edge
%% Initializing Variables 
dx = 0.1; % grid discrete distance
dy = 0.1; % grid discrete distance

Nx = x/dx; % x direction grid
Ny = y/dy; % y direction grid
 
phi = zeros (Nx*Ny,1); % phi initialization
phiOld = zeros (Nx*Ny,1); % phiOld initialization
miu = zeros (Nx*Ny,1); % miu initialization

% PDE Coefficient initialization
a = zeros (Nx*Ny,1);
b = zeros (Nx*Ny,1);
c = zeros (Nx*Ny,1);
d = zeros (Nx*Ny,1);
e = zeros (Nx*Ny,1);
f = zeros (Nx*Ny,1);
g = zeros (Nx*Ny,1);
 
% j = -1; % y direction tracing
% i = -1; % x direction tracing
%  
u_ = -1; % local phi derivative to x initialization
v_ = -1; % local phi derivative to y initialization
 
m = -1; % Mech number locally
A = zeros (Nx*Ny,1); % initialize A matrix
 
error = Inf; % Begining with infinite error
tol = 1E-4; % Tolerence
count = 0; % Runtime counting

xspan = linspace(0,x,Nx); % x discrete spacing
dydx = toc * (-4 * xspan + 82); % Dy/Dx
dydx(xspan<xle | xspan>xte ) = 0; % Zero Dy/Dx outside the airfoil range 

errorlist = nan(1,1000); % Storing error

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
            A(loc,1) = (1 - Minf^2) - (gamma + 1) * Minf^2 ...
                / Uinf * (phi(loc+1,1) - phi(loc-1,1) ) / ( 2 * dx);
            if A(loc,1) > 0
                miu(loc,1) = 0;
            else
                miu(loc,1) = 1;
            end
        end
    end
    % Storing phiOld for error computing
    phiOld = phi;
    % factor
    for j = 2:(Ny-1)
        for i = 3:(Nx-1)
            % Locating specific Phi position
            loc = (j-1) * Ny + i;
            % Coefficient of (i,j-1) 
            c(loc,1) = 1 / dy^2;
            % Coefficient of (i-2,j)
            g(loc,1) = miu(loc-1,1) * A(loc-1,1) / dx^2;
            % Coefficient of (i-1,j)
            d(loc,1) = (1-miu(loc,1)) * A(loc,1) / (dx)^2 ...
                - 2 * miu(loc-1,1) * A(loc-1,1) / (dx)^2;
            % Coefficient of (i+1,j)
            e(loc,1) = (1-miu(loc,1)) * A(loc,1) / (dx)^2;
            % Coefficient of (i,j+1))
            b(loc,1) = 1 / dy^2;
            % Coefficient of (i,j)
            a(loc,1) = miu(loc-1,1) * A(loc-1,1) / (dx)^2 - 2 * ((1-miu(loc,1))*A(loc,1))/(dx)^2 - 2/(dy)^2;
            % Phi Computing
            phi(loc,1) = -(c(loc,1) * phi(loc-Nx,1) ...
                + g(loc,1) * phi(loc-2,1) ...
                + d(loc,1) * phi(loc-1,1) ...
                + e(loc,1) * phi(loc+1,1) ...
                + b(loc,1) * phi(loc+Nx,1) ...
                ) ./ a(loc,1); 
        end
    end

    error = max(max(abs(phiOld - phi)));
    count = count + 1;
    errorlist (count) = error;
end
%% plot
% Preparing Phi into [Nx,Ny] Matrix
plotphi = zeros(Nx,Ny); % Converted Phi initialization
for j = 1 : Ny
    for i = 1 : Nx
        loc = (j-1) * Ny + i;
        plotphi(i,j) = phi(loc,1);
    end
end

% Problem 2
% Error
semilogy(errorlist)
xlabel('Iterations')
ylabel('$L_{\infty}$','interpreter','latex')
str = join({'Convergence history ','( Mach = ',num2str(Minf),' )'});
title(str)


% Cp Computing
cp = zeros (Nx,Ny); % cp initialization
 for i = 2:Nx-1
     for j = 1:Ny
        u_ = (plotphi(i+1,2) - plotphi(i-1,2))/(2*dx) ;
        cp (i,j)  = 2*u_/Uinf;
     end
 end
 xaxis = 0:dx:50;
 figure (21)
 plot (cp);
 xlim([19.5/dx,21.5/dx]);
 ylim([-0.5,1]);
 
 % Pressure Contour
 p = nan(Nx,Ny);
 for i = 2:Nx-1
     for j = 1:Ny
        u_ = (plotphi(j,i+1) - plotphi(j,i-1))/(2*dx) ;
        p(i,j) = Pinf * (1 - gamma*Minf^2 * (u_ /Uinf));
     end 
 end
 figure (22)
 contourf(p)
 colorbar
xlim([19.5/dx,21.5/dx]);
ylim([0,10]);
 set(gca,'PlotBoxAspectRatio',[2 1 1]);

  figure (23)

 contourf(plotphi)

 
  dx = 0.2;
