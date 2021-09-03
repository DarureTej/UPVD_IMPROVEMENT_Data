


     syms    Tcout    Tcin    Toa   Fc   I   Ttout Ttin Ft

          %% solar collector Model selection 
            UL_1  =   3.89;        % solar heat loss coefficient (W m-2 K-1)
            UL_2  =   0.01;        % solar heat loss coefficient (W m-2 K-1)
            eta   =   0.77;        % efficiency (dimentionaless)
            Ac    =   9.38;        % solar collector plate surface area (m2)
            rho   =   1000;        % solar collector fluid density (kg m-3)
            c     =   4180;        % Solar collector plate specific heat (J kg °C)
            V     =  0.003;        % Solar collector fluid volume (m3)


            % Estimation of coefficents 
            C   = rho*c*V;

            a1  = Ac*eta/C;     % coefficent of I
            a2  = -UL_1*Ac/C;    % coefficent of loss
            a3  = -UL_2*Ac/C;     % coefficent of loss
            a4  = V^-1;         % coefficent of FC*Tcin
            a5  = -a4;          % coefficent of FC*Tcout



            % State equation for solar collector (Refer to Notes_Solar_hot_water_system for more details)  

            state_1    =    a1*I + a2*((Tcin+Tcout)/2-Toa) + a3*((Tcin+Tcout)/2-Toa)^2 + a4*Fc*Tcin + a5*Fc*Tcout;
            %% Hot Water Tank Model

            A_t      =      0.5;            % area of hot water tank for heat exchange (m2)
            U_t      =      250;            % solar heat loss coefficient (W m-2 K-1)
            rho      =     1000;            % water density (kg m-3)
            c        =     4200;            % water specific heat (J kg °C)
            V_ct     =    0.075;            % solar fluid vomule to heat transfer (m3)
            V_t      =    0.075;            % hot water tank volume (m3)

            C_t    =  rho*c*V_t;    

            % Defining the variables 

          
            % State equation for hot water tank (Refer to Notes_Solar_hot_water_system for more details)  

            state_2    =  (V_ct^-1)*Fc*(Tcout-Tcin) +(U_t*A_t/2/C_t)*(Tcout-Ttin)-(U_t*A_t/2*C_t)*(Ttout-Tcin);

            state_3    =  (V_t^-1)*Ft*(Ttin-Ttout) +(U_t*A_t/2/C_t)*(Ttin-Tcout)+(U_t*A_t/2*C_t)*(Tcin-Ttout);
            
  
           A_solar_tank=jacobian([state_1;state_2;state_3],[Tcout Tcin Ttout]);

           B_solar_tank=jacobian([state_1;state_2;state_3],[Fc Ft I Toa Ttin]);   
           
           
           %% Operating point values from [J. Buzas, 'Modeling and simulation fo solat thermal system', MAC in Simulation, 1998 ]

             % states OP
            Tcout    =     40; 
            Tcin     =     16.87; 
            Ttout    =     60;    

            % Input OP
            Fc       =     0.15;
            Ft       =     0.15;    

            % dist OP
            Ttin     =     10; 
            Toa      =     30; 
            I        =     1350;                                                    

        %% Model data

            A_solar_tank_1      =     eval(A_solar_tank);

            B_solar_tank_1      =     eval(B_solar_tank(:,1:2));

            Bd_solar_tank_1     =     eval(B_solar_tank(:,3:5));

            C_solar_tank_1      =     eye(3);

            D_solar_tank_1      =     zeros(3,2);