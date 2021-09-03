function dxdt = dynamics_SHWS(t,x,u,d)
 % This function defines the dynamics of the solar hot water system with counter current heat exchanger. 
 % Tcout is output temperature of solar collector fluid (°C)
 % Tcin is input temperature  of solar collector fluid (°C)
 % Ttout is output water temperature of heat exchanger  (°C)
 % Tcin is input water temperature of heat exchanger (°C)
 % Fc is flow rate of solar collector fluid (m3 s-1)
 % Ft is water flow rate in heat exchanger (m3 s-1)
 % Toa is outside tempetaure (°C)
 % I is solar irradiation on solar plate collectors (Wm-2)

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

            C_c      = rho*c*V_c;
            C_t      =  rho*c*V_t;   
            C_ct     =  rho_c*c_c*V_ct;   
            
            
% dxdt = zeros(3,1);          
                       
Tcout = x(1);
Tcin = x(2);
Ttout = x(3);

Fc= u(1);
Ft=u(2);
      
I = d(1);
Toa =d(2);
Ttin= d(3);

% system equations 
dxdt(1,1)    =  (A_c*eta/C_c)*I + (-U_c*A_c/C_c)*((Tcin+Tcout)/2-Toa) + Fc/V_c*(Tcin-Tcout);           

type_heat_exchanger = 2; 

switch type_heat_exchanger
   
    case 0
    % Cocourant
    dxdt(2,1)    =  0;
    dxdt(3,1)    =  0; 
    
    case 1
    % Cocourant
    dxdt(2,1)    =  (V_ct^-1)*Fc*(Tcout-Tcin) -(U_t*A_t/C_ct)*(Tcin-Ttout);
    dxdt(3,1)    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcin); 
        
    case 2
    %Contre-courant
    dxdt(2,1)    =  (V_ct^-1)*Fc*(Tcout-Tcin) -(U_t*A_t/C_ct)*(Tcin-Ttin);
    dxdt(3,1)    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcout); 
        
end


end