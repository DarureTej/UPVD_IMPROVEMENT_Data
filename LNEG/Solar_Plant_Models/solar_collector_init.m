% Solar collector model initialization data 
clear;clc
% Acknoweldgment: This work contains technical information created under Project IMPROVEMENT 
% funded by European Commission with the European Regional Development Funds (ERDF) under the program Interreg SUDOESOE3/P3/E0901.
% Author: Tejaswinee Darure Date: 20 Jan, 2020
%% solar collector specifications 

            U_c      =          7;              % solar heat loss coefficient (W m-2 K-1)
            eta      =        1.2;              % optical efficiency (dimentionaless)
            A_c       =         2;              % solar collector plate surface area (m2)
            rho_c    =       1043;              % solar collector fluid density (kg m-3)
            c_c      =       4180;              % Solar collector plate specific heat (J kg °C)
            V_c      =     0.0075;              % Solar collector fliud volume (m3)

% Estimation of coefficents 
           C_c      = rho_c*c_c*V_c;

           
            % State equation for solar collector 

    
%             state_1    =     (A_c*eta/C_c)*I + (-U_c*A_c/C_c)*((Tcin+Tcout)/2-Toa) + Fc/V_c*(Tcin-Tcout);
            
            
            a1  = A_c*eta/C_c;          % coefficent of I
            a2  = -U_c*A_c/C_c/2;       % coefficent of Tcin and Tcout
            a3  = U_c*A_c/C_c;          % coefficent of Toa
            a4  = V_c^-1;               % coefficent of FC*Tcin
            a5  = -a4;                  % coefficent of FC*Tcout
      
            
%% Heat Exchanger specifications 

            A_t      =      0.5;              % area of hot water tank for heat exchange (m2)
            U_t      =      250;            % solar heat loss coefficient (W m-2 K-1)
            rho      =     1000;            % water density (kg m-3)
            c        =     4200;            % water specific heat (J kg °C)
            V_ct     =    0.075;            % solar fluid volume to heat trannsfer (m3)
            V_t      =    0.075;            % hot water tank volume (m3)

            C_t      =  rho*c*V_t;   
            C_ct     =  rho_c*c_c*V_ct;   
