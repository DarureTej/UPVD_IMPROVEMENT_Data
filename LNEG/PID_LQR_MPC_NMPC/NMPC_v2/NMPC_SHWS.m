clear, clc, close all

Num_days = 30;                       % Number of days: 90 ie. three months (1 particular season)
tsim = 24*Num_days;                  % Per minute for 90 days period
N  = 30;                             % Temperature over a day is sinusoidal
nu = 2;                              % number of inputs
ts = 3600;                             % sampling time 


udata = zeros(2,1);                 % data initialization 
ddata = [1000;20;16];               % data initialization 
x = [60;27.5;47.5];
xdata  = x;
xref_plt = repmat([65;51.5],1,tsim);

u_min = [0;0];                      % lower boundary condition 
u_max = [2.5e-5;2.1e-5];            % lower boundary condition 
u_opt =repmat([1.4e-5;0.9e-5],1,N);

%% system simulation

for time = 1: 1: tsim-N
    disp(' ')
    disp(['time stamp ',num2str(time),'/',num2str(tsim-1)])
    
    
    
    disturbance_case = 1;          % [1 = varying disturbance, 2 = constant disturbance]
    
    switch disturbance_case
        case 1
            [distp,distt]=dist_predicition_N(time,N,2);
        case 2
            distt  = [1000;20;10];
            distp  =  repmat(distt,1,N);
    end
   
    %  NMPC  Controller
    options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxFunEvals',10000,'MaxIter',1000);
    
    [u_opt,fobj,exitflag] = fmincon(@(u) objective_SHWS(x,u,distp,N),[u_opt(:,2:end) u_opt(:,end)],[],[],[],[],repmat(u_min,N,1),repmat(u_max,N,1),[],options);
    
    %  Application of controller output on the system
    u = u_opt(:,1);
    [time,x_plus] = ode45(@(t,x) dynamics_SHWS(t,x,u,distt), [0,ts], x);
    x =  x_plus(end,:)';
    
    % saving data
    xdata  = [xdata x];
    udata = [udata u ];
    ddata= [ddata distt];
    
    
end


Filename = ['Data_',string(Num_days),'Days_MPC',datestr(now,'dd_mmm'),'.mat']
save (Filename,'xdata','udata','ddata','xref_plt(:,1:tsim-N+1)');
%% plotting the results 
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

function [fob] = objective_SHWS(x,u_opt,distp,N)
% [fob] = objective_SHWS(xref,wx,wu)% 
% Input dimentions are given as below:
% xref denoted the references and are explained as xref  = [Tcout-setpoint, Ttout_setoint]
% wx are the weights on the states and are given as xref = [Tcout-weight, Ttout_weight]
% wu are weights on inputs and given as wu = [Tcout-weight, Ttout_weight]

fob = 0;                                 % initialization
ts  = 3600;                              % sampling time 
x_N =[];                                 % initialization
xref = repmat([65 51.5],N,1);            % reference setpoints for x(1), x(3)
wx = [1e0,1e0];                          % weights for states 
wu = [1e1;1e1]*0;                        % weights for inputs
 
 % calculation of objective function
    for i = 2: N
         d= distp(:,i);
         [t,x_pred]  = ode45_dynamics(@dynamics_SHWS,u_opt(i-1:i),d,[0 ts],x);
          x_N = [x_N x_pred(end,:)'];
%           fob = fob + wx(1)*norm((xref(1)-x_pred(end,1)), 2)  + wx(2)*norm((xref(2)-x_pred(end,3)), 2)+   wu(1)*norm(u_opt(i-1), 2)+   wu(2)* norm(u_opt(i), 2);
    end

     fob = wx(1)*norm((xref(:,1) - x_N(1,:)), 2) + wx(2)*norm((xref(:,2) - x_N(3,:)), 2);
end

