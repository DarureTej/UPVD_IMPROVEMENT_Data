function [sys,x0,str,ts] = Systeme(t,x,e,flag,x10,x20,x30)
%
% S-function pour simuler le modèle d'un turbo-réacteur
%
%
% Paramètres passés à la fonction : 
%  x10, x20, x30 : valeurs initiales des états

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(x10,x20,x30);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,e);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,e);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 2, 4, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end
% end csfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes(x10,x20,x30)

sizes = simsizes;
sizes.NumContStates  = 3;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 4;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [x10 x20 x30] ;           % Valeurs initiales de l'état
str = [];
ts  = [0 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,e)

Fc = e(1);
Ft = e(2);
I = e(3);
Toa = e(4);
% Ttin = e(5);
% 
% Tcout = x(1);
% Tcin = x(2);
% Ttout = x(3);

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

%Ancien code
% dTcout    =  (A_c*eta/C_c)*I + (-U_c*A_c/C_c)*((Tcin+Tcout)/2-Toa) + Fc/V_c*(Tcin-Tcout);           
% dTcin    =  (V_ct^-1)*Fc*(Tcout-Tcin) +(U_t*A_t/C_ct)*(Ttout-Tcin);
% dTtout    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcin); 

Ttin = e(5);
Tcout = x(1);
Tcin = x(2);
Ttout = x(3);

%Panneau solaire
dTcout    =  (A_c*eta/C_c)*I + (-U_c*A_c/C_c)*((Tcin+Tcout)/2-Toa) + Fc/V_c*(Tcin-Tcout);           

type_echangeur = 2;

switch type_echangeur
   
    case 0
    % Cocourant
    dTcin    =  0;
    dTtout    =  0; 
    
    case 1
    % Cocourant
    dTcin    =  (V_ct^-1)*Fc*(Tcout-Tcin) -(U_t*A_t/C_ct)*(Tcin-Ttout);
    dTtout    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcin); 
        
    case 2
    %Contre-courant
    dTcin    =  (V_ct^-1)*Fc*(Tcout-Tcin) -(U_t*A_t/C_ct)*(Tcin-Ttin);
    dTtout    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcout); 
        
end

sys = [dTcout dTcin dTtout];

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,e)
% 
% Tcout = x(1);
% Tcin = x(2);
% Ttout = x(3);
%sys = [Tcout Tcin Ttout];

Ttin = e(5);
Tcout = x(1);
Tcin = x(2);
Ttout = x(3);

sys = [Ttin Tcout Tcin Ttout];

% end mdlOutputs