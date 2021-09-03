% Solar collector model initialization data 

% Acknoweldgment: This work contains technical information created under Project IMPROVEMENT 
% funded by European Commission with the European Regional Development Funds (ERDF) under the program Interreg SUDOESOE3/P3/E0901.
% Author: Tejaswinee Darure Date: 20 Jan, 2020
%% solar collector specifications 

UL  =  7;        % solar heat loss coefficient (W m-2 K-1)
eta = 0.8;       % optical efficiency (dimentionaless)
Ac  =  2;        % solar collector plate surface area (m2)
rho = 1043;      % solar collector fluid density (kg m-3)
c   = 4180;      % Solar collector plate specific heat (J kg °C)
V   = 0.075;     % Solar collector fliui volume (m3)


% Estimation of coefficents 
C   = rho*c*V;

a1  = Ac*eta/C;     % coefficent of I
a2  = UL*Ac/2/C;    % coefficent of Tcin and Tcout
a3  = -UL*Ac/C;     % coefficent of Ta
a4  = V^-1;         % coefficent of FC*Tcin
a5  = -a4;          % coefficent of FC*Tcout