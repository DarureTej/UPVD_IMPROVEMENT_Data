clear, clc, close all

Num_days = 1;            % Number of days: 90 ie. three months (1 particular season)
tsim = 90;% 60*24*Num_days;   % Per minute for 90 days period
N  = 30;                 % Temperature over a day is sinusoidal
nu = 2;                  % number of inputs
ts = 60;                 % sampling time 


% Initialization 
udata = zeros(2,1);
ddata = [1000;20;16];
x = [60;27.5;47.5];
xdata  = x;
xref_plt = repmat([65;51.5],1,tsim);
% boundary condition 
u_min = [0;0];
u_max = [2.5e-5;2.1e-5];
u_opt =repmat([1.4e-5;0.9e-5],1,N);
%% system simulation 
for time = 1: 1: tsim-N
    
    disp(' ')
    disp(['time stamp ',num2str(time),'/',num2str(tsim-1)])
    
    
  
    disturbance_case = 2;          % [1 = varying disturbance, 2 = constant disturbance]
    
    switch disturbance_case
        case 1
            [distp,distt]=dist_predicition_N(t,N);
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
%%
subplot(3,1,1)
plot(xdata')

hold on 
plot(xref_plt')
legend('Tcout (°C)', 'Tcin (°C)','Ttout (°C)','Tcout-ref(°C)','Ttout-ref(°C)');grid; 

subplot(3,1,2)
plot(udata')
legend('Fc(m3s-1)', 'Ft(m3s-1)');
grid

subplot(3,1,3)
plot(ddata(1,:)')
legend('I(W/m3)');
grid