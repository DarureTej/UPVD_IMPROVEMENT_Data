%% Mathematical Model of entire system

% Acknowledgment: This work contains technical information created under Project IMPROVEMENT 
% funded by European Commission with the European Regional Development Funds (ERDF) under the program Interreg SUDOESOE3/P3/E0901.
% Author: Tejaswinee Darure Date: 04 Feb, 2021

clear;clc
%% solar collector Model 

U_L1      =          3.89;              % solar heat loss coefficient (W m-2 K-1)
U_L2      =          0.01;              % solar heat loss coefficient (W m-2 K-1)
eta      =        0.77;              % optical efficiency (dimentionaless)
Ac       =          9.38;              % solar collector plate surface area (m2)
rho_c    =       1000;              % solar collector fluid density (kg m-3)
c_c      =       4180;              % Solar collector plate specific heat (J kg °C)
V_c      =     0.003;              % Solar collector fliud volume (m3)


% Estimation of coefficents 

C_c      =     rho_c*c_c*V_c;

a1       =     Ac*eta/C_c;          % coefficent of I
a2       =      -U_L1*Ac/2/C_c;        % coefficent of Tcin and Tcout
a3       =     -U_L2*Ac/C_c;         % coefficent of Ta
a4       =     V_c^-1;              % coefficent of FC*Tcin
a5       =     -a4;                 % coefficent of FC*Tcout


syms Tcout Tcin Toa Fc I 
state_1    = (Ac*eta/C_c)*I + (-U_c*Ac/C_c)*Toa + (U_c*Ac/2/C_c)*Tcout + ((U_c*Ac/2/C_c))*Tcin + (V_c^-1)*Fc*Tcin + (-V_c^-1)*Fc*Tcout;

%% Hot Water Tank Model
A_t   =  0.5;          % area of hot water tank for heat exchange (m2)
U_t   =  7;            % solar heat loss coefficient (W m-2 K-1)
rho   = 1000;          % water density (kg m-3)
c     = 4200;          % water specific heat (J kg °C)
V_ct   = 0.075;        % solar fluid vomule to heat trannsfer (m3)
V_t   = 0.075;         % hot water tank volume (m3)

C_t   =  rho*c*V_t;    

syms Ttout Ttin Ft

state_2    =  (V_t^-1)*Fc*(Tcout-Tcin) +(-U_t*A_t/C_t)*(Tcin-Ttout);
state_3    =  (V_t^-1)*Ft*(Ttin-Ttout) +(-U_t*A_t/C_t)*(Ttout-Tcin);


%% Thermal Storage Tank Model
V_s   = 0.075;         % volumne of thermal storage tank (m3)

C_t   =  rho*c*V_s; 

syms Tsout Tsin Fhe Qhp 

state_4 = (V_s^-1)*Ft*Ttout + Qhp/C_t+Fhe*(Tsout-Tsin)/V_s;


%% Heat Exchanger - FCU 

k       = 0.1  ;            % coefficient of conductive heat transfer for air (dimentionless)
h       = 0.2 ;            % coefficient of convective heat transfer for air (dimentionless)
Ao      =  1   ;       % shell-side area (m2) 
delta_z = 0.05 ;             % length of the tube (dimentionless)
rho_a   = 1000 ;        % air density (kg m-3)
c_a     = 0.5  ;            % heat capacity of air 
V_he    = 0.075;        % FCU-heat exchanger volume (m3)
Ts      =  45;             % surface temperature (degC)

C_he   =  k*Ao/(rho*c*V_he*delta_z);

syms Taout Tain Fa

state_5 = (V_he^-1)*Fa*(Tain-Taout)+(-C_he)*(Taout-Tsin);
state_6 = (V_he^-1)*Fa*(Tsout-Tsin)+(-C_he)*(Tsin-Taout)+ (- h*Ao'/rho/c/V_he)*(Ts - Toa);

%% Calculation of Jacobian matrix
A=jacobian([state_1;state_2;state_3;state_4;state_5;state_6],[Tcout Tcin Ttout Tsout Tsin Taout]);
B=jacobian([state_1;state_2;state_3;state_4;state_5;state_6],[Fc Ft Fhe Fa Ttin Toa I]);
C=jacobian( [Tcout Tcin Ttout Tsout Tsin Taout],[Tcout Tcin Ttout Tsout Tsin Taout]);
D=jacobian([Tcout Tcin Ttout Tsout Tsin Taout],[Fc Ft Fhe Fa Ttin Toa I]);

% A_solar_tank=jacobian([state_1;state_2;state_3],[Tcout Tcin Ttout]);
% B_solar_tank=jacobian([state_1;state_2;state_3],[Fc Ft Fhe]);

% finding operating point 
Tcout=40; Tcin=16.87; Ttout=60; Tsout=40; Tsin=20; Taout=45; Tain =25;      % states OP
Fc=1.5e-15; Ft=1.5e-15; Fhe=1.e-5; Fa=1.e-5;                                % input OP
Ttin=10; Toa=30; I=1350;                                                    % dist OP

% 
% A_solar_tank_1 = eval(A_solar_tank)
% B_solar_tank_1 = eval(B_solar_tank)
A1=eval(A);B1=eval(B);C1=eval(C);D1=eval(D);
% sys1=c2d(ss(A1,B1(:,1:4),C1,D1(:,1:4)),1);
%  
% % Adjusting state space matrices 
% A_bidon = -eye(7,7)+eye(7,7)*0.5+0.03*rand(7,7); cam
E1 = eye(7);
sys = ss(A1,B1,C1,D1,1)

figure(1)
subplot(1,2,1)
spy(A1)
subplot(1,2,2)
spy(B1)
