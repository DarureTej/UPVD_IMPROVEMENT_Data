% Solar collector model initialization data 

% Acknoweldgment: This work contains technical information created under Project IMPROVEMENT 
% funded by European Commission with the European Regional Development Funds (ERDF) under the program Interreg SUDOESOE3/P3/E0901.
% Author: Tejaswinee Darure Date: 20 Jan, 2020

%Updated by Mariana Silva Date: 6 April, 2021
%% solar collector specifications 

UL_1  =  3.89;        % solar heat loss coefficient (W m-2 K-1)
UL_2 = 0.01;        % solar heat loss coefficient (W m-2 K-1)
eta = 0.77;       % efficiency (dimentionaless)
Ac  =  9.38;        % solar collector plate surface area (m2)
rho = 1000;      % solar collector fluid density (kg m-3)
c   = 4180;      % Solar collector plate specific heat (J kg °C)
V   = 0.003;     % Solar collector fliui volume (m3)



% Estimation of coefficents 
C   = rho*c*V;

a1  = Ac*eta/C;     % coefficent of I
a2  = -UL_1*Ac/C;    % coefficent of loss
a3  = -UL_2*Ac/C;     % coefficent of loss
a4  = V^-1;         % coefficent of FC*Tcin
a5  = -a4;          % coefficent of FC*Tcout