%%% MPC program for 6 zones with economic/energy minimization :
clc
clear all
yalmip('clear')
close all 
%    _ _ _ _ _ _ _ _  _ _ _ _ _ _  _ _ _ _ _ _ _  _ _ _ _ _  _ _
%   |           |                          |                     |
%   |           |           5              |                     |
%   |           |                          |                     | 
%   |     4     |_ _ _ _ _ _ _ _ _ _ _ _  _|                     | 
%   |                                      |           6         | 
%   |_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ |                     |
%   |           |           |              |                     |
%   |      1    |    2      |     3        |                     |   
%   |_ _ _ _ _ _|_ _ _ _ _ _| _ _ _ _ _ _  |_ _ _ _ _ _ _ _ _ _ _|  
%                   
%           LNEG Lab Building layout - 5 rooms 
% [1,2,3 --> offices, 4 --> corridor, 5 --> meeting room, 6--> open office space]

%% import the dynamic mode for 6 zones- building 
run('jacobianmodel_6zones_FCU');
A=sys.a ;Bu1=sys.b(:,1:nu);Bu2=sys.b(:,nu+1:2*nu+1); Bd=sys.b(:,2*nu+1:end); C=sys.c; 
tsim=24*4*5;

%% decomposiiton into subsystems in overalapping structure 
A11=[A(1:2,1:2) A(1:2,5);A(5,1:2) A(5,5)]; A12=[A(1:2,3:4);A(5,3:4)]; B1T=[Bu1(1:2,1:2) Bu1(1:2,5);Bu1(5,1:2) Bu1(5,5)];
Bd1=[Bd(1:2,1:2) Bd(1:2,5) Bd(1:2,7);Bd(5,1:2) Bd(5,5) Bd(5,7)];B1F=[Bu2(1:2,1:2) Bu2(1:2,5);Bu2(5,1:2) Bu2(5,5)];
A22=[A(3:4,3:4) A(3:4,6);A(6,3:4) A(6,6)];A21=[A(3:4,1:2);A(6,1:2)];B2T=[Bu1(3:4,3:4) Bu1(3:4,6);Bu1(6,3:4) Bu1(6,6)];
Bd2=[Bd(3:4,3:4) Bd(3:4,6) Bd(3:4,7);Bd(6,3:4) Bd(6,6) Bd(6,7)];B2F=[Bu2(3:4,3:4) Bu2(3:4,6);Bu2(6,3:4) Bu2(6,6)];
xl = sdpvar(1,1);xu = sdpvar(1,1);

run('generate_reference') % generate setpoint-comfort range 
nd=nx+1;
Np =24;  % horizon as one day
uop=ones(nx/2,1)*25;
k1 = [1 1 1]*0.5;  % cost per kw for central supply fan 
k2= 1; % cost per kj by central heating coil
R = diag(ones(1,nu/2)*1); % decide imporatnce of inputs as chnage in airflow will be at priority 

%%  controller for zone 1
% Initial state
uu1 = sdpvar(repmat(nu/2,1,Np),ones(1,Np));
xx1 = sdpvar(repmat(nx/2,1,Np),ones(1,Np));
lambda1 = sdpvar(repmat(3,1,Np-1),ones(1,Np-1));
xx34= sdpvar(2,1);zeta1 = sdpvar(1,1);
dd1 = sdpvar(repmat(4,1,Np),ones(1,Np));
constraints1 = [];objective1 = 0;
uf1=ones(3,1)*0.192;
for k = 1:Np-1
 objective1=objective1+k1*((uu1{k}+uop)).^2+k2*sum(uu1{k}+uop)+zeta1^2;
 constraints1 = [constraints1, xx1{k+1} == A11*xx1{k}+A12*xx34+B1T*uu1{k}+Bd1*dd1{k}+B1F*uf1];
 constraints1 = [constraints1, -3 <=uu1{k}<= 3,xl-zeta1<=xx1{k}<=xu+zeta1,zeta1 >= 0,zeta1 <= 0.5 ];
end
parameters_in1={xx1{1},[dd1{:}],xx34,xl,xu};
solutions_out1 = {[uu1{:}], [xx1{:}],objective1};
options = sdpsettings('verbose',1);
controller1 = optimizer(constraints1, objective1,options,parameters_in1,solutions_out1);

%%  controller for zone 2 
% Initial state
uu2 = sdpvar(repmat(nu/2,1,Np),ones(1,Np));
xx2 = sdpvar(repmat(nx/2,1,Np),ones(1,Np));
dd2 = sdpvar(repmat(4,1,Np),ones(1,Np));
xx12= sdpvar(2,1);zeta2 = sdpvar(1,1);
constraints2 = [];objective2 = 0;
uf2=ones(3,1)*0.192;
for k = 1:Np-1
 objective2 = objective2 + k1*((uu2{k}+uop)).^2+k2*sum(uu2{k}+uop)+zeta2^2;
constraints2 = [constraints2, xx2{k+1} == A22*xx2{k}+A21*xx12+B2T*uu2{k}+Bd2*dd2{k}+B2F*uf1];
constraints2 = [constraints2, -3 <= uu2{k}<= 3,xl-zeta2<=xx2{k}<=xu+zeta2,zeta2 >= 0,zeta2 <= 0.5];
end
parameters_in2={xx2{1},[dd2{:}],xx12,xl,xu};
solutions_out2 = {[uu2{:}], [xx2{:}],objective2};

controller2 = optimizer(constraints2, objective2,options,parameters_in2,solutions_out2);
enerd=zeros(1,tsim);
%% Implemenatation on building 
xx=zeros(nx,1);xx1=zeros(nx/2,1);xx2=zeros(nx/2,1);xx12=zeros(2,1);xx34=zeros(2,1);
for t=2:tsim 
% % load disturbance 
    [distp,distt,distds]=distpredicition6zone(t,Np);      
%  
% % solving an optimization problem for subsystem 1
    [solutions1,diagnostics1] = controller1{{xx1,distds,xx34,Ulb(t),Uub(t)}};
     if diagnostics1 == 1
        error('The problem is infeasible');
        break;
     end
     U1 = solutions1{1}; X1 = solutions1{2};     
%      
%  % solving an optimization problem for subsystem 1
    [solutions2,diagnostics2] = controller2{{xx2,distds,xx12,Ulb(t),Uub(t)}};
     if diagnostics2 == 1
        error('The problem is infeasible');
        break;
     end
     U2 = solutions2{1}; X2 = solutions2{2};  
   
    
% % Applying on building 
  xx=sys.a*xx+sys.b(:,1:2*nu)*[U1(1:2,1);U2(1:2,1);U1(3,1);U2(3,1);uf1;uf2]+Bd*distt; 
  y=C*x+0.007*rand(nx,1);
  
% % save data
% 
UU(:,t)=[U1(1:2,1);U2(1:2,1);U1(3,1);U2(3,1)]+[uop;uop];
XX(:,t)=[xx+ones(nx,1)*23];
YY(:,t)=y+ones(nx,1)*23;
Dist(:,t)= distt+[ones(nx,1)*0.65 ;Toa];
enerd(:,t)=enerd(:,t-1)+sum((UU(:,t)).^2)+sum(UU(:,t));
xx1=[xx(1:2);xx(5)];
xx2=[xx(3:4);xx(6)];
xx12=xx(1:2);xx34=xx(3:4);
end 
xl=Ulb(1:tsim)+22.5;xu=Uub(1:tsim)+23.5;

%% title('Temperature for every zone')
figure(1)
subplot(6,1,1)
plot((1:tsim)/25/4,XX(1,:)','Color',[0.2,0.2,1],'lineWidth',3)
ylim([15 30]);grid;title('Temperature for Office 1');xlabel('Days');ylabel('Temp(°C)');
hold on;plot((1:tsim)/25/4,xl,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');hold on;plot((1:tsim)/25/4,xu,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');
subplot(6,1,2)
plot((1:tsim)/25/4,XX(2,:)','Color',[0,0.2,1],'lineWidth',3)
ylim([15 30]);grid;title('Temperature for Office 2');xlabel('Days');ylabel('Temp(°C)');
hold on;plot((1:tsim)/25/4,xl,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');hold on;plot((1:tsim)/25/4,xu,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');
subplot(6,1,3)
plot((1:tsim)/25/4,XX(3,:)','Color',[0,0.2,1],'lineWidth',3)
ylim([15 30]);grid;title('Temperature for Office 3');xlabel('Days');ylabel('Temp(°C)');
hold on;plot((1:tsim)/25/4,xl,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');hold on;plot((1:tsim)/25/4,xu,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');
subplot(6,1,4)
plot((1:tsim)/25/4,XX(4,:)','Color',[0,0.2,1],'lineWidth',3)
ylim([15 30]);grid;title('Temperature for Corridor');xlabel('Days');ylabel('Temp(°C)');
hold on;plot((1:tsim)/25/4,xl,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');hold on;plot((1:tsim)/25/4,xu,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');
subplot(6,1,5)
plot((1:tsim)/25/4,XX(5,:)','Color',[0,0.2,1],'lineWidth',3)
ylim([15 30]);grid;title('Temperature for Meeting room ');xlabel('Days');ylabel('Temp(°C)');
hold on;plot((1:tsim)/25/4,xl,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');hold on;plot((1:tsim)/25/4,xu,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');
subplot(6,1,6)
plot((1:tsim)/25/4,XX(6,:)','Color',[0,0.2,1],'lineWidth',3)
ylim([15 30]);grid;title('Temperature for Open Office Space');xlabel('Days');ylabel('Temp(°C)');
hold on;plot((1:tsim)/25/4,xl,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');hold on;plot((1:tsim)/25/4,xu,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');
%% title('Control inputs')
figure(2)
subplot(6,1,1)
plot((1:tsim)/25/4,UU(1,:)','Color',[0,0.5,0.8],'lineWidth',3)
ylim([0 40]);grid;title('Supply air Tempetaure rate for Office 1');xlabel('Days');ylabel('Temp(°C)');
subplot(6,1,2)
plot((1:tsim)/25/4,UU(2,:)','Color',[0,0.5,0.8],'lineWidth',3)
ylim([0 40]);grid;title('Supply air Tempetaure  for Office 1');xlabel('Days');ylabel('Temp(°C)');
subplot(6,1,3)
plot((1:tsim)/25/4,UU(3,:)','Color',[0,0.5,0.8],'lineWidth',3)
ylim([0 40]);grid;title('Supply air Tempetaure  for Office 1');xlabel('Days');ylabel('Temp(°C)');
subplot(6,1,4)
plot((1:tsim)/25/4,UU(4,:)','Color',[0,0.5,0.8],'lineWidth',3)
ylim([0 40]);grid;title('Supply air Tempetaure  for Corridor');xlabel('Days');ylabel('Temp(°C)');
subplot(6,1,5)
plot((1:tsim)/25/4,UU(5,:)','Color',[0,0.5,0.8],'lineWidth',3)
ylim([0 40]);grid;title('Supply air Tempetaure  for Meeting room');xlabel('Days');ylabel('Temp(°C)');
subplot(6,1,6)
plot((1:tsim)/25/4,UU(6,:)','Color',[0,0.5,0.8],'lineWidth',3)
ylim([0 40]);grid;title('Supply air Tempetaure  for  Open Office Space');xlabel('Days');ylabel('Temp(°C)');
%%
figure(3)
subplot(2,1,1)
plot((1:tsim)/25/4,Dist(1,:)','Color',[0,0.5,0],'lineWidth',3)  
grid;title('Heat flux due to occupancy ');xlabel('Days');ylabel('heat flux (kW)');
subplot(2,1,2)
plot((1:tsim)/25/4,Dist(7,:)'+2.5,'Color',[0,0,0.5],'lineWidth',3.5)
grid;xlabel('Days');ylabel('Temperature(°C)');title('Weather Temperature ');
%

%%
figure(4)
subplot(2,1,1)
plot((1:tsim)/25/4,XX(1,:)','Color',[0,0.5,0],'lineWidth',3)  
ylim([15 30]);grid;title('Temperature in Office 1 ');xlabel('Days');ylabel('Temperature(°C) ');
hold on;plot((1:tsim)/25/4,xl,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');hold on;plot((1:tsim)/25/4,xu,'Color',[1,0.5,0],'lineWidth',2,'LineStyle','--');
subplot(2,1,2)
plot((1:tsim)/25/4,UU(1,:)','Color',[0,0,0.5],'lineWidth',3.5)
grid;xlabel('Days');ylabel('Temperature(°C)');title('Supply Air Temperature for Office 1 ');
ylim([19 26])


