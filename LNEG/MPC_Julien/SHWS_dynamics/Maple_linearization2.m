clc;clear;

%% defining all the variables and parameters
syms U_c eta A_c rho_c c_c V_c
syms A_t U_t rho c  V_ct V_t
syms  Tcout    Tcin    Toa   Fc   I
syms Ttout Ttin Ft

C_c   = rho*c*V_c;
C_t   =  rho*c*V_t;
C_ct  =  rho_c*c_c*V_ct;

%% differential equations for solar collector and heat exchager

state_1    =  (A_c*eta/C_c)*I + (-U_c*A_c/C_c)*((Tcin+Tcout)/2-Toa) + Fc/V_c*(Tcin-Tcout);

type_echangeur = 2;

switch type_echangeur
    
    case 0
        % Cocourant
        dTcin    =  0;
        dTtout    =  0;
        
    case 1
        % Cocourant
        state_2    =  (V_ct^-1)*Fc*(Tcout-Tcin) -(U_t*A_t/C_ct)*(Tcin-Ttout);
        state_3    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcin);
        
    case 2
        %Contre-courant
        state_2    =  (V_ct^-1)*Fc*(Tcout-Tcin) -(U_t*A_t/C_ct)*(Tcin-Ttin);
        state_3    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcout);
        
end

states = [state_1;state_2;state_3];

% Equ = solve(states,[Tcin,Fc,Ft]);

A_solar_tank=jacobian([state_1;state_2;state_3],[Tcout Tcin Ttout])

B_solar_tank=jacobian([state_1;state_2;state_3],[Fc Ft I Toa Ttin])

%%
% syms x1(t) x2(t) x3(t)
% ode1 = diff(x1)    ==  (A_c*eta/C_c)*I + (-U_c*A_c/C_c)*Toa + (U_c*A_c/2/C_c)*Tcout + ...
%                                        ((U_c*A_c/2/C_c))*Tcin + (V_c^-1)*Fc*Tcin + (-V_c^-1)*Fc*Tcout;
%
% ode2 = diff(x2)    ==  (V_ct^-1)*Fc*(Tcout-Tcin) +(-U_t*A_t/C_ct)*(Tcin-Ttout);
%
% ode3 =diff(x3)    ==  (V_t^-1)*Ft*(Ttin-Ttout) +(-U_t*A_t/C_t)*(Tcin-Ttout);
%
% odes = [ode1;ode2;ode3];
%
% cond1 = x1(0) == 0;
% cond2 = x2(0) == 1;
% cond3 = x3(0) == 1;
% conds = [cond1; cond2; cond2];
%
% S = dsolve(odes,conds);


%% substituting values


U_c      =          7;              % solar heat loss coefficient (W m-2 K-1)
eta      =        1.2;              % optical efficiency (dimentionaless)
A_c       =         2;              % solar collector plate surface area (m2)
rho_c    =       1043;              % solar collector fluid density (kg m-3)
c_c      =       4180;              % Solar collector plate specific heat (J kg °C)
V_c      =     0.0075;              % Solar collector fliud volume (m3)

A_t      =      0.5;              % area of hot water tank for heat exchange (m2)
U_t      =      250;            % solar heat loss coefficient (W m-2 K-1)
rho      =     1000;            % water density (kg m-3)
c        =     4200;            % water specific heat (J kg °C)
V_ct     =    0.075;            % solar fluid volume to heat trannsfer (m3)
V_t      =    0.075;            % hot water tank volume (m3)

C_c   = rho*c*V_c;
C_t    =  rho*c*V_t;
C_ct    =  rho_c*c_c*V_ct;

%% OP

switch type_echangeur
    case  1 % cocourrant
        Tcout    =     89;
        Tcin     =     63.5;
        Ttout    =     50;
        
        % Input OP
        Fc       =     1.5e-5;
        Ft       =     1.0e-5  ;
        
        % dist OP
        Ttin     =     10;
        Toa      =     20;
        I        =     1000;
        
        
    case  2 % contre courant 
        Tcout    =     60;
        Tcin     =     47.5;
        Ttout    =     27.5;
        
        % Input OP
        Fc       =     1.5e-5;
        Ft       =     1.0e-5  ;
        
        % dist OP
        Ttin     =     10;
        Toa      =     20;
        I        =     1000;
        
end
%% Model data

A_solar_tank_1      =     eval(A_solar_tank);

B_solar_tank_1      =     eval(B_solar_tank(:,1:2));

Bd_solar_tank_1     =     eval(B_solar_tank(:,3:5));

C_solar_tank_1      =     eye(3);

D_solar_tank_1      =     zeros(3,5);

%% controllability check
nx      =      3;                             % Number of states
nu      =      2;                             % Number of inputs
nd      =      3;                             % Number of disturbances


eig(A_solar_tank_1)
if eig(A_solar_tank_1)<0
    display('------- System is stable------- ')
elseif eig(A_solar_tank_1)<=0
    display('------- System is marginally stable------- ')
else
    display('------- System is UNSTABLE------- ')
end


if nx == rank(ctrb(A_solar_tank_1,B_solar_tank_1))
    display('------- System (A, B) is controllable------- ')
else
    display('------- System (A, B) is NOT controllable -------')
end

if nx == rank(ctrb(A_solar_tank_1,B_solar_tank_1(:,1)))
    display('------- System (A, B(:,1)) is controllable considering only first (Fc) input------ ')
else
    display('------- System (A, B) is NOT controllable considering only first (Fc) input------ ')
end


if nx == rank(ctrb(A_solar_tank_1,B_solar_tank_1(:,2)))
    display('------- System (A, B(:,2)) is controllable considering only second (Ft) input------ ')
else
    display('------- System (A, B) is NOT controllable considering only second (Ft) input------ ------- ')
end

if nx == rank(obsv(A_solar_tank_1,C_solar_tank_1))
    display('------- System (A, C) is observable------- ')
else
    display('------- System (A, B) is NOT observable -------')
end


% sys = ss(A_solar_tank_1,[B_solar_tank_1 Bd_solar_tank_1],C_solar_tank_1,D_solar_tank_1);
% display('Bd= ')
% sys.b(:,3:end)
% step(sys)