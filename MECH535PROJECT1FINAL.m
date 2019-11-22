%%%                 MECH 309 - Numerical Methods in Mech Eng                 %%%

% Presented to Prof Siva Nadarajah Winter 2019 - November 22th

%Yiming Yao 260769906
%Zechen Ren 260765431
%

clear all
close all
clc
for COMP=1:-1:0 %This will provide compressible solution followed by incompressible
% COMP = 0; %SET TO ZERO FOR INCOMPRESSIBLE/ZERO LOSS DESIGN

%%                 GIVEN PARAMETERS                 %%

% Before rotor, rC?1 = 39.3m2/s (constant with radius: free vortex)
% After rotor, rC?2 = 117.8m2/s (constant with radius: free vortex)
% N = 6000rpm
% rshroud = 0.50m
% rhub = 0.45m
% mass flow(m) = 30.4kg/s
% To(inlet) = 288K
% po(inlet) = 1.5kg/m3
% w, at blade trailing edge = 0.03

%%                 SETTING KNOWN VARIBLES                 %%

OMEGA = 6000*pi*2/60; %[rad/s]
RSHROUD = 0.5; %[m]
RHUB = 0.45; %[m]
XMASS = 30.4; %[kg/s]
TO1 = 288; %[K]
PO1 = 100000; %[Pa]
DENSO1 = 1.5; %[kg/m3]
if COMP==1
OMEGLOS = 0.3; %[@ T.E.]
else
OMEGLOS = 0;    
end  
CP = 1005; %[J/kg*K]
GAMMA = 1.4/(1.4-1);
RGAS = 287; %[J/kg*K]
NSTATN = 51;
NSTRM = 11;
NLE=21;
NTE=31;
RCUBR=39.3;  %m2/s
RCUAR=117.8; %m2/s
H = (RSHROUD-RHUB)/(NSTRM-1);
delr=(RSHROUD-RHUB)/(NSTRM-1);
delz=delr;

%%                 Initializing Variables with Dimension                   %%
                        %In the following form: 
%   [1,1  1,2  1,3  1,4  1,5 ....  1,51
%    2,1  2,2  2,3  2,4  2,5 ....  2,51
%    3,1  3,2  3,3  3,4  3,5 ....  3,51
%      .    .    .    .    . ....     .
%   11,1 11,2 11,3 11,4 11,5 .... 11,51]

PSI = zeros(NSTRM,NSTATN);
DENSITY = ones(NSTRM,NSTATN);
RADIUS = zeros(NSTRM,NSTATN);
RCU = zeros(NSTRM,NSTATN);
CZ = zeros(NSTRM,NSTATN);
CR = zeros(NSTRM,NSTATN);

HTOTAL = zeros(NSTRM,NSTATN);
PTOTAL = zeros(NSTRM,NSTATN);

RHS = zeros(NSTRM,NSTATN);
ENTROPY = zeros(NSTRM,NSTATN);

A = zeros(NSTRM,NSTATN);
B = zeros(NSTRM,NSTATN);
DENRAD = zeros(NSTRM,NSTATN);

%%                 Setting Tolerance                 %%
TOLPSI = 0.00001;
TOLRHS = 0.01;
TOLDEN = 0.001;

%%                 1. INITIALIZE VALUE OF PSI, DENSITY, RADIUS, HTOTAL, PTOTAL, RCU                 %%

%Inputing initial HTOTAL values 
HTOTAL(:,:) = TO1*CP;
PTOTAL(:,:) = PO1;

for j=0:10
    RADIUS(j+1,:)= RSHROUD-j*delr;
end
for j=1:NSTRM
    PSI(j,:) = (((RADIUS(j,j)).^2)-(RHUB.^2))/((RSHROUD.^2)-(RHUB.^2));
end
DENSITY(:,:) = DENSO1;

% initializing RCU which is given in the project and unchanging   
for i=1:NSTATN
    if (i <= NLE)
        RCU(:,i)=RCUBR; 
    elseif (i >= NTE)                   
        RCU(:,i)=RCUAR;   
    else
        RCU(:,i)=RCUBR+((RCUAR-RCUBR)/(NTE-NLE))*(i-NLE);   
    end       
end

%%                 2. Initializing CR & CX                 %%

%CR stays the zero as with the inital values of psi there is zero gradient 
%in the axial direction as such CR is zero
for i=1:NSTATN
for j=NSTRM:-1:1
    if j==NSTRM
        CZ(j,i)=(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j-1,i)-PSI(j,i))/delr);
    elseif j==1
        CZ(j,i)=(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j,i)-PSI(j+1,i))/delr);
    else
        CZ(j,i)=(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j-1,i)-PSI(j+1,i))/(2*delr));
    end
end
end

%%                 6. REPEAT STEP 2-5 UNTIL CONVERGENCE (Big Loop)                 %%
                            %%% START OF MAIN LOOP %%%
ERRDEN=1; %To get loop started
ERRRHS=1;

NN=0; %COUNTER FOR MAIN LOOP

while ((ERRRHS>TOLRHS && ERRDEN>TOLDEN) && NN<100)
      NN=NN+1;

    %%                 2. CALCULATING THE VELOCITY FIELD CX AND CR                 %%

    for i=1:NSTATN
    for j=NSTRM:-1:1
        if i==1
           CR(j,i)=-(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j,i+1)-PSI(j,i))/delz);
        elseif i==NSTATN
           CR(j,i)=-(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j,i)-PSI(j,i-1))/delz);
        else
           CR(j,i)=-(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j,i+1)-PSI(j,i-1))/(2*delz));
        end
    end
    end   
    for i=1:NSTATN
    for j=NSTRM:-1:1
        if j==NSTRM
            CZ(j,i)=(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j-1,i)-PSI(j,i))/delr);
        elseif j==1
            CZ(j,i)=(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j,i)-PSI(j+1,i))/delr);
        else
            CZ(j,i)=(XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i)))*((PSI(j-1,i)-PSI(j+1,i))/(2*delr));
        end
    end
    end

    %%                 3. CALCULATING THERMODYNAMIC VARIABLES STATION BY STATION                 %%
        
    % 3.1 UPDATING THE DENSITY ON THE INLET
        for j=1:NSTRM
            CSQ = ((RCU(j,1)/RADIUS(j,1)).^2) + CZ(j,1).^2 + CR(j,1).^2;
            HSTATIC = HTOTAL(j,1) - CSQ/2;
            PSTATIC = PTOTAL(j,1) * (HSTATIC/HTOTAL(j,1)).^GAMMA;
            if COMP == 1
            DENSITY(j,1) = PSTATIC/(RGAS*HSTATIC/CP);
            end
        end

    % 3.2 SWEEP THE PLANES FROM PLANE 2 TO NSTTN

        NLE =21;
        NTE =31;

        ERRDENS = 0;
        CHANGE = 0;

        for i=2:NSTATN
            
                LEFT = i-1;
            if (i > NLE) && (i < NTE)
                LEFT = NLE;
            end

            NSTART = NSTRM;
            PSIDN = PSI(NSTART,LEFT);
            PSIUP = PSI(NSTART-1,LEFT);

            % Locating Origins Streamline %
            for j=NSTRM:-1:1
               DESIRED = PSI(j,i);
                while (~(DESIRED <= PSIUP && DESIRED >= PSIDN))
                    NSTART = NSTART-1;
                    PSIDN = PSI(NSTART,LEFT);
                    PSIUP = PSI(NSTART-1,LEFT);
                    if ( NSTART == 1 )
                        disp('CANNOT TRACE STREAMLINE')
                        break
                    end
                end
                % streamline origin is found 
                DELTA = (DESIRED-PSIDN)/(PSIUP-PSIDN);

                ROTATE = 0;
                if (i > NLE) && (i <= NTE)
                    ROTATE = OMEGA;
                end

                RCU1       =DELTA*(RCU(NSTART-1,LEFT)-RCU(NSTART,LEFT))+RCU(NSTART,LEFT);
                HTOTAL1    =DELTA*(HTOTAL(NSTART-1,LEFT)-HTOTAL(NSTART,LEFT))+HTOTAL(NSTART,LEFT);
                PTOTAL1    =DELTA*(PTOTAL(NSTART-1,LEFT)-PTOTAL(NSTART,LEFT))+PTOTAL(NSTART,LEFT);
                RAD1       =DELTA*(RADIUS(NSTART-1,LEFT)-RADIUS(NSTART,LEFT))+RADIUS(NSTART,LEFT);

                DENS1      =DELTA*(DENSITY(NSTART-1,LEFT)-DENSITY(NSTART,LEFT))+DENSITY(NSTART,LEFT);
                CZ1        =DELTA*(CZ(NSTART-1,LEFT)-CZ(NSTART,LEFT))+CZ(NSTART,LEFT);
                CR1        =DELTA*(CR(NSTART-1,LEFT)-CR(NSTART,LEFT))+CR(NSTART,LEFT);
                ENTROP1    =DELTA*(ENTROPY(NSTART-1,LEFT)-ENTROPY(NSTART,LEFT))+ENTROPY(NSTART,LEFT);
                C1SQ       =CZ1.^2+CR1.^2+(RCU1/RAD1).^2;
                HSTAT1     =HTOTAL1 - (C1SQ/2);
                PSTAT1     =PTOTAL1 * (HSTAT1/HTOTAL1).^GAMMA ;

                ROTALP1    =HTOTAL1 - ROTATE * RCU1;
                HOR2       =ROTALP1 + ((ROTATE * RADIUS(j,i)).^2 )/ 2;
                HOR1       =ROTALP1 + (ROTATE * RAD1).^2 / 2;
                POR1       =PTOTAL1 * (HOR1/HTOTAL1).^GAMMA;
                POR2IDL    =POR1 * (HOR2/HOR1).^GAMMA;

                if (i > NLE) && (i <= NTE)
                    OMEGLOSS = OMEGLOS*(i-NLE)/(NTE-NLE);
                    PLOSS = OMEGLOSS * (POR1-PSTAT1);
                else
                    OMEGLOSS = 0.0;
                    PLOSS = 0.0;
                end
                POR2       =POR2IDL - PLOSS;

                %PROJECT 1: ROTOR WITH RCU SPECIFIED
                RCU(j,i) = RCU1;           
                if ( i > NLE && i <= NTE )
                    RCU(:,i)=RCUBR+((RCUAR-RCUBR)/(NTE-NLE))*(i-NLE);   
                end
                
                %COMMON CALCULATION BLOCK

                CU          =RCU(j,i)/RADIUS(j,i);
                VU          =CU - ROTATE * RADIUS(j,i);
                V2SQ        =VU.^2 + CZ(j,i).^2 + CR(j,i).^2;
                C2SQ        =CU.^2 + CZ(j,i).^2 + CR(j,i).^2;
                HSTATIC     =HOR2-V2SQ/2;
                HTOTAL(j,i) =HSTATIC + C2SQ/2;
                PTOTAL(j,i) =POR2 * (HTOTAL(j,i)/HOR2).^GAMMA;
                PSTATIC     =PTOTAL(j,i) * (HSTATIC/HTOTAL(j,i)).^GAMMA;
                if COMP == 1
                DENSOLD     =DENSITY(j,i);
                DENSITY(j,i)=PSTATIC/(RGAS*HSTATIC/CP);
                CHANGE      =abs(DENSOLD-DENSITY(j,i));
                ERRDENS     =max(CHANGE,ERRDENS);
                end
                ENTROPY(j,i)=CP * log(HTOTAL(j,i)/HTOTAL1) - RGAS * log(PTOTAL(j,i)/PTOTAL1) + ENTROP1;
            end
        end

    %%                 4. ASSEMBLE THE RHS FOR UNKNOWN NODES I.E. ALL NODES EXCEPT HUB/SHROUD                 %%

        ERRRHS = 0.0;
        CHANGE = 0.0;
        for i=1:NSTATN
            for j=NSTRM-1:-1:2
                RHSOLD      =RHS(j,i);
                CU          =RCU(j,i)/RADIUS(j,i);
                TSTATIC     =(HTOTAL(j,i)-(CZ(j,i).^2+CR(j,i).^2+CU.^2)/2)/CP;
                RHS(j,i)    = -(1/CZ(j,i)) * (2*pi/(2*delr*XMASS)) * (CU/RADIUS(j,i)*...
                              (RCU(j-1,i)-RCU(j+1,i))-(HTOTAL(j-1,i)-HTOTAL(j+1,i))+...
                              (TSTATIC*(ENTROPY(j-1,i)-ENTROPY(j+1,i))));
                CHANGE      =abs(RHS(j,i)-RHSOLD);
                ERRRHS      =max(ERRRHS,CHANGE);
            end
        end

    %%                 5. UPDATE THE STREAM FUNCTION PSI TO CONVERGENCE                 %%

        ERRPSI=1; %SET TO START LOOP
        N=1;
        DENRAD = DENSITY.*RADIUS;

        while ERRPSI>TOLPSI && N<500
                N=N+1;
                ERRPSI=0;
                for i=1:NSTATN
                    for j=NSTRM-1:-1:2
                        if i == NSTATN
                            B(j,i) = (2*PSI(j,i-1)/(DENRAD(j,i)))+...
                                     (PSI(j-1,i)/((DENRAD(j-1,i)+DENRAD(j,i))/2))+...
                                     (PSI(j+1,i)/((DENRAD(j+1,i)+DENRAD(j,i))/2));
                            A(j,i) = 1/((2/(DENRAD(j,i)))+...
                                     (1/((DENRAD(j-1,i)+DENRAD(j,i))/2))+...
                                     (1/((DENRAD(j+1,i)+DENRAD(j,i))/2)));                 
                        elseif i == 1
                            B(j,i) = (2*PSI(j,i+1)/(DENRAD(j,i)))+...
                                     (PSI(j-1,i)/((DENRAD(j-1,i)+DENRAD(j,i))/2))+...
                                     (PSI(j+1,i)/((DENRAD(j+1,i)+DENRAD(j,i))/2));
                            A(j,i) = 1/((2/(DENRAD(j,i)))+...
                                     (1/((DENRAD(j-1,i)+DENRAD(j,i))/2))+...
                                     (1/((DENRAD(j+1,i)+DENRAD(j,i))/2)));                              
                        else
                            B(j,i) = (PSI(j,i+1)/((DENRAD(j,i+1)+DENRAD(j,i))/2))+...
                                     (PSI(j,i-1)/((DENRAD(j,i-1)+DENRAD(j,i))/2))+...
                                     (PSI(j-1,i)/((DENRAD(j-1,i)+DENRAD(j,i))/2))+...
                                     (PSI(j+1,i)/((DENRAD(j+1,i)+DENRAD(j,i))/2));
                            A(j,i) = 1/((1/((DENRAD(j,i+1)+DENRAD(j,i))/2))+...
                                     (1/((DENRAD(j,i-1)+DENRAD(j,i))/2))+...
                                     (1/((DENRAD(j-1,i)+DENRAD(j,i))/2))+...
                                     (1/((DENRAD(j+1,i)+DENRAD(j,i))/2)));
                        end
                        PSIOLD=PSI(j,i);
                        PSIPSI = A(j,i)*(B(j,i)+(delz.^2)*RHS(j,i));
                        PSI(j,i)=0.5*PSIPSI+(1-0.5)*PSI(j,i);
                        CHANGE = abs(PSI(j,i)-PSIOLD);
                        ERRPSI = max(CHANGE,ERRPSI);     
                    end
                end
        end
end

%%                  OUTPUT OF RESULTS: QUESTION 1 & 2                 %%
k=0;
for i=1:NSTATN
    Z = real(i-1)*H;
    k=k+1;
    disp('STATION:'); disp(k)
    disp('  RADIUS:      PSI:    CX:        CR:      CM:       CU:       VU:       C:        V:')
    for j=1:NSTRM
        CU = RCU(j,i)/RADIUS(j,i);
        ROTATE = 0;
        if ( i >= NLE && i <= NTE); ROTATE = OMEGA;end
        VU = CU-ROTATE*RADIUS(j,i);
        C = sqrt(CR(j,i).^2+CZ(j,i).^2+CU.^2);
        V = sqrt(CR(j,i).^2+CZ(j,i).^2+VU.^2);
        CM = sqrt(CR(j,i).^2+CZ(j,i).^2);

        disp([RADIUS(j,i),PSI(j,i),CZ(j,i),CR(j,i),CM,CU,VU,C,V])
    end
    disp('  PSTATIC:  TSTATIC:  TTOTAL:   TREL:     PTOTAL:   PREL:      BETAZ:    ALPHAZ:    ')
    for j=1:NSTRM
        CU = RCU(j,i) / RADIUS(j,i);
        ROTATE = 0;
        if ( i >= NLE && i <= NTE); ROTATE = OMEGA;end
        VU = CU - ROTATE*RADIUS(j,i);
        BETAZ = atan(VU/CZ(j,i))*(180/pi);
        if i == NTE
            BETA_TE (j,1) = BETAZ;
        end
        ALPHAZ = atan(CU/CZ(j,i))*(180/pi);
        CSQ = CR(j,i).^2+CZ(j,i).^2+CU.^2;
        VSQ = CR(j,i).^2+CZ(j,i).^2+VU.^2;
        TTOTAL = HTOTAL(j,i)/CP;
        TSTATIC = TTOTAL - CSQ/(2*CP);
        TREL = TSTATIC + VSQ/(2*CP);
        PSTATIC = DENSITY(j,i)*RGAS*TSTATIC;
        PREL = PTOTAL(j,i)*(TREL/TTOTAL).^GAMMA;
        
        disp([PSTATIC/1000,TSTATIC,TTOTAL,TREL,PTOTAL(j,i)/1000,PREL/1000,BETAZ,ALPHAZ])    
    end
    disp('  RCU:       USPEED:   DENSITY:   Z:        ENTROPY:')
    for j=1:NSTRM
       ROTATE = 0;
       if ( i >= NLE && i <= NTE); ROTATE = OMEGA;end
       USPEED = RADIUS(j,i)*ROTATE;
       
       disp([RCU(j,i),USPEED,DENSITY(j,i),Z,ENTROPY(j,i)]) 
    end    
end
%%                   PLOTS: QUESTION 3                      %%
if COMP == 1
figure(1)
plot (RADIUS(:,NTE),CZ(:,NTE))
xlabel('Radius [m]');
ylabel('Axial velocity [m/s]');
title('Compressible flow with losses: Axial velocity vs. Radius');

figure(2)
plot (RADIUS(:,NTE),CR(:,NTE))
xlabel('Radius [m]');
ylabel('Radial velocity [m/s]');
title('Compressible flow with losses: Radial velocity vs. Radius');


figure(3)
plot (RADIUS(:,NTE),BETA_TE)
xlabel('Radius [m]');
ylabel('Blade Angle [degree]');
title('Compressible flow with losses: Blade angle vs. Radius');

figure(4)
plot (RADIUS(:,NTE),DENSITY(:,NTE))
xlabel('Radius [m]');
ylabel('Density [kg/s]');
title('Compressible flow with losses: Density vs. Radius');
else
figure(5)
plot (RADIUS(:,NTE),CZ(:,NTE))
hold on
CZ_ANA = 136*ones(NSTRM,1);
plot (RADIUS(:,NTE),CZ_ANA)
xlabel('Radius [m]');
ylabel('Axial velocity [m/s]');
title('Incompressible flow without losses: Axial velocity vs. Radius');
legend('Numerical','Analytical');

figure(6)
plot (RADIUS(:,NTE),CR(:,NTE),'LineWidth',4)
hold on
CR_ANA = zeros(NSTRM,1);
plot (RADIUS(:,NTE),CR_ANA,'LineWidth',2)
xlabel('Radius [m]');
ylabel('Radial velocity [m/s]');
title('Incompressible flow without losses: Radial velocity vs. Radius');
legend('Numerical','Analytical');

figure(7)
plot (RADIUS(:,NTE),BETA_TE,'LineWidth',4)
hold on
BETA_ANA = -427.2*RADIUS(:,NTE)+183.6;
plot (RADIUS(:,NTE),BETA_ANA,'LineWidth',2)
xlabel('Radius [m]');
ylabel('Blade Angle [degree]');
title('Incompressible flow without losses: Blade angle vs. Radius');
legend('Numerical','Analytical');

figure(8)
plot (RADIUS(:,NTE),DENSITY(:,NTE),'LineWidth',4)
hold on
DENSITY_ANA = 1.5*ones(NSTRM,1);
plot (RADIUS(:,NTE),DENSITY_ANA,'LineWidth',2)
xlabel('Radius [m]');
ylabel('Density [kg/s]');
title('Incompressible flow without losses: Density vs. Radius');
legend('Numerical','Analytical');
end

%%                 3D BLADE SHAPE: QUESTION 4                 %%

BETADEG=zeros(NTE-NLE+1,NTE-NLE+1);
YMAG=zeros(NTE-NLE+1,NTE-NLE+1);
BLADESHAPE=zeros(NTE-NLE+1,NTE-NLE+1);

for i=NLE:NTE
    for j=1:NSTRM
        CU = RCU(j,i) / RADIUS(j,i);
        ROTATE = 0;
        if ( i >= NLE && i <= NTE); ROTATE = OMEGA;end
        VU = CU - ROTATE*RADIUS(j,i);
        BETADEG(j,i-NLE+1) = atan(VU/CZ(j,i));
        YMAG(j,i-NLE+2)=CZ(j,i)*tan(BETADEG(j,i-NLE+1))*abs(VU/CZ(j,i))*delz; %abs(VU/CZ(j,i))*delz is the scaling factor
    end
end

for i=2:11
BLADESHAPE(:,i)=YMAG(:,i)+BLADESHAPE(:,i-1);
end

%Data from analytical solution @LE beta1h=55.15 beta2s=60  @TE betalh=8.64 beta2s=30

BETA_ANA=zeros(11,11);
BETA_ANA_LE = linspace(60,55.15,11)';
BETA_ANA_TE = linspace(30,8.64,11)';

for j=1:11
    BETA_ANA(j,:)=-linspace(BETA_ANA_LE(j),BETA_ANA_TE(j),11);
end

BETA_ANA(:,:) = BETA_ANA(:,:)*pi/180;

 for i=1:11
    for j=1:11   
        CU = RCU(j,i+20) / RADIUS(j,i);
        ROTATE=OMEGA;
        VU = CU - ROTATE*RADIUS(j,i);
        YANA(j,i+1)=136*tan(BETA_ANA(j,i))*abs(VU/136)*delz;
    end
end

BLADEANA = zeros(11,11);

for i=2:11
BLADEANA(:,i)=YANA(:,i)+BLADEANA(:,i-1);
end




if COMP == 1
figure(9)
legend('Compressible')
s=surf(BLADESHAPE);
s.EdgeColor = 'none';
s.FaceColor = 'r';
title('Blade Shape with BetaZ')
xlabel('Leading Edge   to   Trailing Edge')
ylabel('RADIUS')
else
hold on
t=surf(BLADESHAPE);
t.EdgeColor = 'none';
t.FaceColor = 'y';
title('Blade Shape with BetaZ')
xlabel('Leading Edge   to   Trailing Edge')
ylabel('RADIUS')
u=surf(BLADEANA);
u.EdgeColor = 'none';
u.FaceColor = 'b';
legend('Incompressible','ANALYTICAL')
hold off
end
end










