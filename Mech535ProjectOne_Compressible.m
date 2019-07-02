%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Project One: Free-Vortex Design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------COMPRESSIBLE CASE--------------------------------
%Kristian Stafford Smith----------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%DUCT%%%%%%%%DUCT%%LE%%%ROTOR%%TE%%%%%DUCT%%%%%%%DUCT%%%%%%%%%%%
%%% ==> 1.........10.........21/////////30..........40.........50  ==>
%%% ==> 1.........10.........21/////////30..........40.........50  ==>
%%% IN  1.........10.........21/////////30..........40.........50  OUT
%%% ==> 1.........10.........21/////////30..........40.........50  ==>
%%% ==> 1.........10.........21/////////30..........40.........50  ==>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%% Initializing Values

RHUB = 0.45; %m
RSHROUD = 0.5; %m
NSTATN=50;%vertical grid lines
NSTRM=10;%number of horizontal grid lines
NLE=21;
NTE=30;
DELTAR=(0.5-0.45)/NSTRM; %m                                           
DELTAZ=DELTAR; %m %assuming r/z=1    
XMASS=30.4; %kg/s
OMEGA=6000*2*pi()/60; %rad/s                                        
GAMMA=1.4;                                       
RGAS=287.15;
CP=1005;
T0=288; %K
DENSITY0=1.5; %kg/m3
P0=DENSITY0*RGAS*T0; %Pa
MAXITS=200;
TOLDENS=0.001;
TOLRHS=.01;
TOLPSI=0.00001;

%Dimensioning Variables
CZ=zeros(NSTRM,NSTATN);
CR=zeros(NSTRM,NSTATN);
RADIUS=zeros(NSTRM,NSTATN);
PSI=zeros(NSTRM,NSTATN);
PSIOLD=zeros(NSTRM,NSTATN);  %%%%needs to reinitialize later for iterations
DENSITY=zeros(NSTRM,NSTATN);
DENSITYOLD=zeros(NSTRM,NSTATN);
RHS=zeros(NSTRM,NSTATN);
RCU=zeros(NSTRM,NSTATN);
HTOTAL=zeros(NSTRM,NSTATN);
PTOTAL=zeros(NSTRM,NSTATN);
ENTROPY=zeros(NSTRM,NSTATN);
OMEGLOS=zeros(NSTRM,NSTATN);
ROTHALPY=zeros(NSTRM,NSTATN);
A=zeros(NSTRM,NSTATN);
B=zeros(NSTRM,NSTATN);
RAD=zeros(NSTRM,NSTATN);
DENS=zeros(NSTRM,NSTATN);
HSTATIC=zeros(NSTRM,NSTATN);
PSTATIC=zeros(NSTRM,NSTATN);
POR=zeros(NSTRM,NSTATN);
HOR=zeros(NSTRM,NSTATN);
BETA=zeros(NSTRM,NSTATN);
PSICHANGE=zeros(NSTRM,NSTATN);

% Initializing RADIUS, PSI, DENSITY, PRESSURE, ENTHALPY-----checked
for i=1:NSTATN %horizontal stations
    for j=1:NSTRM %vertical stations
        if j==1                                       
            RADIUS(j,i)=0.45;
        else
            RADIUS(j,i)=0.45+0.05*j/NSTRM; %metres 
        end
            PSI(j,i)=((RADIUS(j,i)^2)-RHUB^2)/(RSHROUD^2-RHUB^2);
            DENSITY(j,i)=1.500; %kg/m3
    end
end

%%%%%%%%%%Set inlet conditions
for j=1:NSTRM
    
  HTOTAL(j,1)=CP*T0;
  PTOTAL(j,1)=P0;
  TTOTAL(j,1)=288;
  ENTROPY(j,1)=0;                   
  
end

% Set and Interpolate Swirl, Loss Coefficient  ----checked2
for j=1:NSTRM
    for i=1:NLE
        RCU(j,i)=39.3; %m2/s
    end
    for i=NLE:NTE
        
        RCU(j,i)=39.3+(117.8-39.3)*(i-NLE)/(NTE-NLE); %m2/s
        
    end
    for i=NTE:NSTATN
        RCU(j,i)=117.8; %m2/s
    end
end


DENSOLD(j,i)=0;
RHSOLD(j,i)=1;

%%%%%%%%%%%%step 3

for w=1:MAXITS
    ERRRHS=0;
   ERRDENS=0;
    

%%%%%%%%%CZ and CR

for i=1:NSTATN
    for j=1:NSTRM
        %%%%for the hub boundary (not corners)
        
       if (i~=1)&&(i~=50)&&(j==1)
           CZ(j,i)=2*XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i))*(PSI(j+1,i)-PSI(j,i))/(2*DELTAR);
            CR(j,i)=2*(-1)*XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i+1)-PSI(j,i-1))/(2*DELTAZ);
        
            %%%shroud boundary
       elseif (i~=1)&&(i~=50)&&(j==10)
               CZ(j,i)=2*XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i)-PSI(j-1,i))/(2*DELTAR);
               CR(j,i)=2*(-1)*XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i+1)-PSI(j,i-1))/(2*DELTAZ);
               
               %%%% at the inlet------should have 2* XMASS?
       elseif (j~=1)&&(j~=10)&&(i==1)
                   CZ(j,i)=XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i))*(PSI(j+1,i)-PSI(j-1,i))/(2*DELTAR);
                    CR(j,i)=(-1)*XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i+1)-PSI(j,i))/(2*DELTAZ);
                    
                    %%%%at outlet boundary
                    
       elseif (j~=1)&&(j~=10)&&(i==50)
                     CZ(j,i)=XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i))*(PSI(j+1,i)-PSI(j-1,i))/(2*DELTAR);  
                    CR(j,i)=(-1)*XMASS/(2*pi*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i)-PSI(j,i-1))/(2*DELTAZ);

%%%%%%%%%%%%% four corners of grid

  %hub left corner
                   elseif (i==1)&&(j==1)
          CZ(j,i) = 2* XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j+1,i)-PSI(j,i))/(2*DELTAR);
          CR(j,i) = 2* (-1)*XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i+1)-PSI(j,i))/(2*DELTAZ);

        %hub right corner
                   elseif	(i==50)&&(j==1)
          CZ(j,i) =2*  XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j+1,i)-PSI(j,i))/(2*DELTAR);
          CR(j,i) = 2* (-1)*XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i)-PSI(j,i-1))/(2*DELTAZ);
 
        %shroud left corner
        elseif(i==1)&&(j==10) 
          CZ(j,i) =2*  XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i)-PSI(j-1,i))/(2*DELTAR);
          CR(j,i) = 2* (-1)*XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i+1)-PSI(j,i))/(2*DELTAZ);
        
        %shroud right corner        
        elseif(i==50)&&(j==10)
          CZ(j,i) =2*  XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i)-PSI(j-1,i))/(2*DELTAR);
          CR(j,i) = 2* (-1)*XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i)-PSI(j,i-1))/(2*DELTAZ);

 % now setting all interior nodes (everywhere else)
      
        else 
          CZ(j,i) = XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j+1,i)-PSI(j-1,i))/(2*DELTAR);
          CR(j,i) = (-1)*XMASS/(2*3.14*DENSITY(j,i)*RADIUS(j,i))*(PSI(j,i+1)-PSI(j,i-1))/(2*DELTAZ);
                   end    
        end
    end

    %% Tracing Thermo Vars
    
    % Update Density of inlet---------checked
    for j=1:NSTRM
        CSQ=(RCU(j,1)/RADIUS(j,1))^2+CZ(j,1)^2+(CR(j,1))^2;
        HSTATIC(j,1)=HTOTAL(j,1)-CSQ/2;
        PSTATIC=PTOTAL(j,1)*(((HSTATIC(j,1)/HTOTAL(j,1)))^(GAMMA/(GAMMA-1)));
        DENSITY(j,1)=PSTATIC/(RGAS*HSTATIC(j,1)/CP);    %for compressible flow, part A! comment out for part b
        
    end

    %Sweep planes
    for i=2:NSTATN
         for j=1:NSTRM
          LEFT=i-1;
             for NSTART=1:9
                PSIDN=PSI(NSTART,LEFT);
                PSIUP=PSI(NSTART+1,LEFT);
                 
                if(PSI(j,i)<=PSIUP)&&(PSI(j,i)>=PSIDN)
                DELTA = (PSI(j,i)-PSIDN)/(PSIUP-PSIDN);

        % LEFT is leading edge for rotor 
        if (i>=21)&&(i<=30)              
            ROTATE=OMEGA;
        else
            ROTATE=0;
        end
        
    %quantities at "1" needed for loss calculation
    RCU1 = DELTA*(RCU(NSTART+1,LEFT)-RCU(NSTART,LEFT))+RCU(NSTART,LEFT);
                HTOTAL1 = DELTA*(HTOTAL(NSTART+1,LEFT)-HTOTAL(NSTART,LEFT))+HTOTAL(NSTART,LEFT);
                PTOTAL1 = DELTA*(PTOTAL(NSTART+1,LEFT)-PTOTAL(NSTART,LEFT))+PTOTAL(NSTART,LEFT);
                RAD1 = DELTA*(RADIUS(NSTART+1,LEFT)-RADIUS(NSTART,LEFT))+RADIUS(NSTART,LEFT);

                %Properties at "1" needed for Loss Calculation
                CZ1 = DELTA*(CZ(NSTART+1,LEFT)-CZ(NSTART,LEFT))+CZ(NSTART,LEFT);
                CR1 = DELTA*(CR(NSTART+1,LEFT)-CR(NSTART,LEFT))+CR(NSTART,LEFT);
                ENTROP1 = DELTA*(ENTROPY(NSTART+1,LEFT)-ENTROPY(NSTART,LEFT))+ENTROPY(NSTART,LEFT);
                C1SQ = CZ1^2+CR1^2+(RCU1/RAD1)^2;
                HSTAT1 = HTOTAL1-C1SQ/2;
                PSTAT1 = PTOTAL1*(HSTAT1/HTOTAL1)^(GAMMA/(GAMMA-1));           %GAMMA

                %Rotating and non-rotating quantities at reference station
                ROTALP1 = HTOTAL1-ROTATE*RCU1;
                HOR2(j,i) = ROTALP1+(1/2)*(ROTATE*RADIUS(j,i))^2;
                HOR1 = ROTALP1+(1/2)*(ROTATE*RAD1)^2;
                POR1 = PTOTAL1*(HOR1/HTOTAL1)^(GAMMA/(GAMMA-1));
                POR2IDL = POR1*(HOR2(j,i)/HOR1)^(GAMMA/(GAMMA-1));
     
    %%%%%%%%%Blade Loss calc --------only for part 1 compressible
       
                    if (i>=NLE) && (i<=NTE)  
                        OMEGLOS = 0.03 * (i-NLE)/(NTE-NLE);
                        PLOSS = OMEGLOS * (POR1 - PSTAT1);
                    else 
                        PLOSS = 0;
                    end
                    
             POR2 = POR2IDL - PLOSS;
        
%%%%%%%%%%%%%%%% COMMON CALC BLOCK

 CU(j,i) = RCU(j,i)/RADIUS(j,i);
                VU = -(ROTATE*RADIUS(j,i) - CU(j,i));
                V2SQ(j,i) = VU^2+CZ(j,i)^2+CR(j,i)^2;
                ALPHA(j,i) = atan (CU(j,i)./sqrt(CZ(j,i)^2+CR(j,i)^2));
                BETA (j,i) = atan (VU/sqrt(CZ(j,i)^2+CR(j,i)^2));
                C2SQ(j,i) = CU(j,i)^2+CZ(j,i)^2+CR(j,i)^2;
                HSTATIC(j,i) = HOR2(j,i)-V2SQ(j,i)/2;            
                HTOTAL(j,i) = HSTATIC(j,i)+C2SQ(j,i)/2;
                TTOTAL(j,i) = HTOTAL(j,i)/CP;
                PTOTAL(j,i) = POR2*(HTOTAL(j,i)/HOR2(j,i))^(GAMMA/(GAMMA-1)); 
                PSTATIC = PTOTAL(j,i)*(HSTATIC(j,i)/HTOTAL(j,i))^(GAMMA/(GAMMA-1));
   %%%%% for compressible
                ERRDENS = max(max(abs(DENSOLD(j,i)-DENSITY(j,i))));
                    DENSOLD (1:NSTRM,1:NSTATN) =  DENSITY(1:NSTRM,1:NSTATN);
                    DENSITY(j,i) = PSTATIC/(RGAS*HSTATIC(j,i)/CP);
    %%%%%%%%
     TSTATIC(j,i)= HSTATIC(j,i)/CP;
     ENTROPY(j,i) = CP*log(HOR2(j,i)/HOR1)-RGAS*log(POR2/POR1)+ENTROP1;
        end
            end
        end
    end

    %%%%%%%%%%%%%%%%%RHS CALC & ERR
    for i = 2:NSTATN 
        for j = 2:(NSTRM-1)
    RHS (j,i) = (2.*pi./(2.*DELTAR* XMASS.*CZ(j,i))).*...
                ((CU(j,i)/RADIUS(j,i)).*((RCU(j+1,i))-(RCU(j-1,i)))+...
                 TSTATIC(j,i).*(ENTROPY(j+1,i)-ENTROPY(j-1,i))-...
                 (HTOTAL(j+1,i)- HTOTAL(j-1,i)));
             
    ERRRHS = max(max(abs(RHSOLD(j,i)-RHS(j,i))));
    RHSOLD (1:NSTRM,1:NSTATN) =  RHS(1:NSTRM,1:NSTATN);
        end
    end       

  %Calc A, B and PSI
          for k = 1: MAXITS
            ERROR = 0;
                for i = 2:NSTATN-1
                    for j = 2:NSTRM-1
                        if i ==  NSTATN-1
                                B(j,i)=   (4*PSI(j,NSTATN-1)./DENSITY(j,i)*RADIUS(j,i)...
                +2*PSI(j+1,i)./((DENSITY(j+1,i)*RADIUS(j+1,i))+(DENSITY(j,i)*RADIUS(j,i)))...
                +2*PSI(j-1,i)./((DENSITY(j-1,i)*RADIUS(j-1,i))+(DENSITY(j,i)*RADIUS(j,i))));
                end
               if j == 2 
               B(j,i)=   (2*PSI(j,i+1)./((DENSITY(j,i+1)*RADIUS(j,i+1))+(DENSITY(j,i)*RADIUS(j,i)))...
                    +2*PSI(j,i-1)./((DENSITY(j,i-1)*RADIUS(j,i-1))+(DENSITY(j,i)*RADIUS(j,i)))...
                    +2*PSI(j+1,i)./((DENSITY(j+1,i)*RADIUS(j+1,i))+(DENSITY(j,i)*RADIUS(j,i))));
               end
               if j == NSTRM-1
               B(j,i)=   (2*PSI(j,i+1)./((DENSITY(j,i+1)*RADIUS(j,i+1))+(DENSITY(j,i)*RADIUS(j,i)))...
                    +2*PSI(j,i-1)./((DENSITY(j,i-1)*RADIUS(j,i-1))+(DENSITY(j,i)*RADIUS(j,i)))...
                    +2*PSI(j+1,i)./((DENSITY(j+1,i)*RADIUS(j+1,i))+(DENSITY(j,i)*RADIUS(j,i)))+1);
               end

                A(j,i)=1./(2./((DENSITY(j,i+1)*RADIUS(j,i+1))+(DENSITY(j,i)*RADIUS(j,i)))...
                          +2./((DENSITY(j,i-1)*RADIUS(j,i-1))+(DENSITY(j,i)*RADIUS(j,i)))...
                          +2./((DENSITY(j+1,i)*RADIUS(j+1,i))+(DENSITY(j,i)*RADIUS(j,i)))...
                          +2./((DENSITY(j-1,i)*RADIUS(j-1,i))+(DENSITY(j,i)*RADIUS(j,i))));
                B(j,i)=   (2*PSI(j,i+1)./((DENSITY(j,i+1)*RADIUS(j,i+1))+(DENSITY(j,i)*RADIUS(j,i)))...
                          +2*PSI(j,i-1)./((DENSITY(j,i-1)*RADIUS(j,i-1))+(DENSITY(j,i)*RADIUS(j,i)))...
                          +2*PSI(j+1,i)./((DENSITY(j+1,i)*RADIUS(j+1,i))+(DENSITY(j,i)*RADIUS(j,i)))...
                          +2*PSI(j-1,i)./((DENSITY(j-1,i)*RADIUS(j-1,i))+(DENSITY(j,i)*RADIUS(j,i))));

                 PSIOLD(j,i) = A(j,i)*(B(j,i)+(DELTAZ^2).*(RHS (j,i)));
                 PSICHANGE(j,i) = abs(PSIOLD(j,i)-PSI(j,i));
                 ERROR = max(max(PSICHANGE));
           end
        end

    PSI(2:NSTRM-1,2:NSTATN-1)=PSIOLD(2:NSTRM-1,2:NSTATN-1);
        
        if ERROR < TOLPSI % ---- Exit condition -----------------
            break;
        end
    end

MAXERR = 0.00001;
TOLDENS=0.001;
TOLRHS=.01;
TOLPSI=0.00001;
% --------- EXIT CONDITION -----------------------------------------
    if ERRRHS < TOLRHS && ERRDENS < TOLDENS
        break;
    end
end

PSI;

%% Plots
%1-D
figure
plot (RADIUS(:,NTE), CZ(:,NTE),RADIUS(:,NTE), CR(:,NTE))
xlabel('RADIUS [m]')
ylabel('Velocity [m/s]')
legend ('CZ', 'CR')
axis([0.45 0.5 -10 170])
grid on
title( 'Velocity vs Radius [COMPRESSIBLE]')

figure
plot (RADIUS(:,NTE), BETA(:,NTE))
xlabel('RADIUS [m]')
ylabel('BETA [radians]')
grid on
title(' Beta vs Radius [COMPRESSIBLE]')

figure
plot (RADIUS(:,NTE), DENSITY(:,NTE))
xlabel('RADIUS [m]')
ylabel('DENSITY [kg/m_3]')
grid on
title('Density vs Radius [COMPRESSIBLE]')

% 3-D

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1,...
    'ZTickLabel',{'Hub','','','','','','','','','Shroud'},...
    'ZTick',[1 2 3 4 5 6 7 8 9 10],...
    'XTickLabel',{'LE','','','','','','','','TE'},...
    'XTick',[1 2 3 4 5 6 7 8 9])

xlabel('Longitudinal Direction');
zlabel('Radial Direction');
box(axes1,'on');
hold(axes1,'all');

BLADE = zeros(NSTRM,(NTE-NLE)); 
for i=1:(NTE-NLE)  
    for j=1:NSTRM
        BLADE(j,i)= BETA(j,NLE+i-1)*DELTAZ;
    end
end

[X,Y]=meshgrid(1:9,1:10);
C=gradient(BLADE);
surf(X,BLADE,Y,C)

        title({'3D Rotor Blade'})

colorbar;       
hold off



