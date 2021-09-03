clear, clc, close all


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
           
            syms    Tcout  Tcin  Toa   Fc   I  Ttout Ttin Ft

            % State equation for solar collector 

            state_1    =     (A_c*eta/C_c)*I + (-U_c*A_c/C_c)*((Tcin+Tcout)/2-Toa) + Fc/V_c*(Tcin-Tcout);
            
           
           % State equation for hot water tank          
            
            type_heat_exchanger = 2; 

            switch type_heat_exchanger

                case 0
                % Cocourant
                state_2    =  0;
                state_3    =  0; 

                case 1
                % Cocourant
                state_2    =  (V_ct^-1)*Fc*(Tcout-Tcin) -(U_t*A_t/C_ct)*(Tcin-Ttout);
                state_3    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcin); 

                case 2
                %Contre-courant
                state_2    =  (V_ct^-1)*Fc*(Tcout-Tcin) -(U_t*A_t/C_ct)*(Tcin-Ttin);
                state_3    =  (V_t^-1)*Ft*(Ttin-Ttout) -(U_t*A_t/C_t)*(Ttout-Tcout); 

            end
           
    A_solar_tank=jacobian([state_1;state_2;state_3],[Tcout Tcin Ttout]);

    B_solar_tank=jacobian([state_1;state_2;state_3],[Fc Ft I Toa Ttin]);
    
    % linearized system at every instant 
    % states OP
    Tcout    =     60; 
    Tcin     =     27.5; 
    Ttout    =     47.5;    

    % Input OP
    Fc       =     1.5e-5;
    Ft       =     1.0e-5;    

    % dist OP
    Ttin     =     10; 
    Toa      =     20; 
    I        =     1000;                                                    

     % state space model 
    A      =     eval(A_solar_tank);
    B      =     eval(B_solar_tank(:,1:2));
    Bd     =     eval(B_solar_tank(:,3:5));
    C      =     [1,0,0;0,0,1];
    D      =     zeros(2,5);
    sysd   =     ss(A,[B Bd],C,D,1);


A = [   -0.0022    0.0018         0;
         0.0002   -0.0006         0;
         0.0004         0   -0.0005];
     
B =    1.0e+03 *[-4.3333         0;
                  0.4333         0;
                        0   -0.5000];
                    
Bd =    1.0e-03 *[0.0762    0.4444         0;
                       0         0    0.3823;  
                       0         0    0.1333];
                   
                   
C =    [ 1     0     0;    0     0     1]; 

D = zeros(2,5);

sysd   =     ss(A,[B Bd],C,D,3600);

Num_days   = 2;                          % Number of days: 90 ie. three months (1 particular season)
tsim       = 24*Num_days;             % Per minute for 90 days period
N          = 10;                         % Temperature over a day is sinusoidal
nu         = 2;                          % number of inputs
ts         = 3600;                          % sampling time 



x_op       = [60;27.5;47.5];             % state operating point 
d_op       = [1000;20;10];               % disturbance initialization
u_op       = [1.5e-5;1.0e-5];            % input initialization

x          = zeros(3,1);                 % state initialization
xref_plt   = repmat([62;49.5],1,tsim);


u_min      = -u_op;                      % lower boundary condition 
u_max      = u_op*2;                     % lower boundary condition 
u_opt      = repmat(zeros(2,1),1,N);

xdata      = x_op;                       % data initialization
udata      = u_op;                       % data initialization 
ddata      = [1000;20;16];               % data initialization 
%% system simulation

for time = 1: 1: tsim-N
    disp(' ')
    disp(['time stamp ',num2str(time),'/',num2str(tsim-1)])
    
    
    
    disturbance_case = 1;                % [1 = varying disturbance, 2 = constant disturbance]
    
    switch disturbance_case
        case 1
            [distp,distt]=dist_predicition_N(time,N,1);
        case 2
            distt  = [1000;20;10];
            distp  =  repmat(distt,1,N);
    end
   
    options   =   optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxFunEvals',10000,'MaxIter',1000);
                                          % NMPC  Controller optimization setting 
    
    [u_opt,fobj,exitflag] = fmincon(@(u) objective_SHWS_MPC(sysd,x,u,distp,N),[u_opt(:,2:end) u_opt(:,end)],[],[],[],[],repmat(u_min,N,1),repmat(u_max,N,1),[],options);
                                          % NMPC  Controller optimization problem 
   
    u          =    u_opt(:,1);                        %  Application of controller output on the system
    [t,x_plus] =    ode45(@(t,x) linearized_fixed_model(t,sysd,x,u,distt), [0 ts], x);
    x          =    x_plus(end,:)';
    
   xdata  = [xdata x+x_op];
   udata = [udata u+u_op ];
   ddata= [ddata distt+d_op];
    
    
end

Filename = ['Data_',string(Num_days),'Days_MPC',datestr(now,'dd_mmm'),'.mat']
save (Filename,'xdata','udata','ddata','xref_plt(:,1:tsim-N+1)');
%% plotting the results 
formatIn = 'dd-mm-yyyy';   
startDate = datenum('03-01-2019',formatIn);
endDate = datenum('01-01-2019',formatIn);
xData = linspace(startDate,endDate,length(xdata));
%% 

figure(1)
subplot(3,1,1)
plot(xData,xdata','lineWidth',1.5)
if Num_days > 30 
datetick('x','mmm-yy','keepticks')
else 
    datetick('x','dd-mmm','keepticks')
end 
hold on 
plot(xData,xref_plt(:,1:tsim-N+1)','lineWidth',1.5);
if Num_days > 30 
datetick('x','mmm-yy','keepticks')
else 
    datetick('x','dd-mmm','keepticks')
end 
legend('Tcout (°C)', 'Tcin (°C)','Ttout (°C)','Tcout-ref(°C)','Ttout-ref(°C)');grid; 

subplot(3,1,2)
plot(xData,udata','lineWidth',1.5)
if Num_days > 30 
datetick('x','mmm-yy','keepticks')
else 
    datetick('x','dd-mmm','keepticks')
end 
legend('Fc(m3s-1)', 'Ft(m3s-1)');
grid

subplot(3,1,3)
plot(xData,ddata(1,:)','lineWidth',1.5)
if Num_days > 30 
datetick('x','mmm-yy','keepticks')
else 
    datetick('x','dd-mmm','keepticks')
end 
legend('I(W/m3)');
grid
%%
function x_plus = linearized_fixed_model(t,sysd,x,u,distt)
   
     A = sysd.A; B = sysd.B(:,1:2); Bd = sysd.B(:,3:5);
    
     x_plus =   A*x+B*u+Bd*distt;
    
end 

function [fob] = objective_SHWS_MPC(sysd,x,u_opt,distp,N)
% [fob] = objective_SHWS(xref,wx,wu)% 
% Input dimentions are given as below:
% xref denoted the references and are explained as xref  = [Tcout-setpoint, Ttout_setoint]
% wx are the weights on the states and are given as xref = [Tcout-weight, Ttout_weight]
% wu are weights on inputs and given as wu = [Tcout-weight, Ttout_weight]

fob = 0; 
ts  = 3600; 
x_N =[];
 % reference setpoints for x(1), x(3)
 xref = repmat([2 0.5],N,1);
 wx = [1e-1,1e-1];
 wu = [1e1;1e1]*0;
 
 % calculation of objective function
    for i = 2: N
         d= distp(:,i);
         u = u_opt(:,1);
         [time,x_pred] = ode45(@(t,x) linearized_fixed_model(t,sysd,x,u,d), [0 ts], x);
         x_N = [x_N x_pred(end,:)'];
%        fob = fob + wx(1)*norm((xref(1)-x_pred(end,1)), 2)  + wx(2)*norm((xref(2)-x_pred(end,3)), 2)+   wu(1)*norm(u_opt(i-1), 2)+   wu(2)* norm(u_opt(i), 2);
    end

     fob = wx(1)*norm((xref(:,1) - x_N(1,:)), 2) + wx(2)*norm((xref(:,2) - x_N(3,:)), 2);
end



