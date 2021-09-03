
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
% [1,2,3 --> offices, 4 --> corridor, 5 --> meeting room, 6--> open office space] _ _ _  _ | 
               

clc;clear;
%%  variables declarartion
nx=6;nu=6;nh=6;
x=sym('x', [nx 1]);    % states == temeprature of each zone
u=sym('u', [nu 1]);    % input == supply air temperature 
lambda=sym('lambda', [nh,1]); % lagrange multiplier 
q=sym('q', [nx 1]); % hea flux due to occupancy 
fsa=sym('fsa', [nu 1]);  
syms Toa % supply air temperature 
%% building properties 
rhoaCpa=1.25625*1.005; % density of air * heat capacity of air
Uw1Aw1=5*9/1000; % heat trnasfer rate for walls shared by both the rooms (1/Rij from equation 2, systol paper)
Uw2Aw2=5*9/1000; % heat trnasfer rate for walls shared by both the rooms (1/Rij from equation 2, systol paper)
URAR=12/1000; % heat trnasfer rate for walls shared by both the roof (1/Rext from equation 2, systol paper)
Cz=47.100; % heat capacity of zone/room (Ci in equation 2 )

h1=(fsa(1)*rhoaCpa*((u(1)-x(1)))+(2*Uw1Aw1*(x(2)-x(1)))+(2*Uw1Aw1*(x(5)-x(1)))+(2*Uw2Aw2*(x(3)-x(1)))+q(1)+(URAR*(Toa-x(1))))/Cz; %energy balance equation for room 1
h2=(fsa(2)*rhoaCpa*((u(2)-x(2)))+(2*Uw1Aw1*(x(1)-x(2)))+(2*Uw1Aw1*(x(5)-x(2)))+(2*Uw2Aw2*(x(4)-x(2)))+q(2)+(URAR*(Toa-x(2))))/Cz;%energy balance equation for room 2
h3=(fsa(3)*rhoaCpa*((u(3)-x(3)))+(2*Uw1Aw1*(x(1)-x(3)))+(2*Uw1Aw1*(x(6)-x(3)))+(2*Uw2Aw2*(x(4)-x(3)))+q(3)+(URAR*(Toa-x(3))))/Cz; %energy balance equation for room 3
h4=(fsa(4)*rhoaCpa*((u(4)-x(4)))+(2*Uw1Aw1*(x(3)-x(4)))+(2*Uw1Aw1*(x(6)-x(4)))+(2*Uw2Aw2*(x(2)-x(4)))+q(4)+(URAR*(Toa-x(4))))/Cz;%energy balance equation for room 4
h5=(fsa(5)*rhoaCpa*((u(5)-x(5)))+(4*Uw1Aw1*(x(2)-x(5)))+(4*Uw2Aw2*(x(1)-x(5)))+q(5)+(URAR*(Toa-x(5))))/Cz; %energy balance equation for room 5
h6=(fsa(6)*rhoaCpa*((u(6)-x(6)))+(4*Uw1Aw1*(x(3)-x(6)))+(4*Uw2Aw2*(x(4)-x(6)))+q(6)+(URAR*(Toa-x(6))))/Cz; %energy balance equation for room 6
h=[h1;h2;h3;h4;h5;h6];


%% Calculation of Jacobian matrix
A=jacobian(h,x);
B=jacobian(h,[u ;fsa ;q  ;Toa ]);
C=jacobian(x,x);
D=jacobian(x,[u ;fsa ;q  ;Toa ]);


%% finding operating point 
Tsa=26; % supply air temperature Ts in equation 2
Toa=5; % outside/weather tempearture 
x1=23;x2=23;x3=23;x4=23;x5=23;x6=23; % states OP initial temperature of romm is 23 deg 
fsa1=0.192;fsa2=0.192;fsa3=0.192;fsa4=0.192;fsa5=0.192;fsa6=0.192; % input OP
u1=25;u2=25;u3=25;u4=25;u5=25;u6=25; 
q1=0.65;q2=0.65;q3=0.65;q4=0.65;q5=0.65;q6=0.65; % dist OP
A1=eval(A);
B1=eval(B);
C1=eval(C);
D1=eval(D);
sys=c2d(ss(A1,B1,C1,D1),60);% sampling time 60sec for discretization
step(sys)
