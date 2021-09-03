load ulisboa.mat
V=10*10*4
ACH=0.3
rho=1.225
CP_air=1.005
Tconf=20
Text=q(2,:);
q4=V*rho*CP_air*ACH*(Text-Tconf)/3.6; % kWh....
q5=q4/2; % In the middle point is half of the inertia....
q(4,:)=q4;
q(5,:)=q5;
q(6,:)=zeros(1,8760);

internal_gains=[0 0 0 0 0 0 0 1000 1000 1000 500 1000 500 500 1000 1000 1000 0 0 0 0 0 0 0]

save ulisboa_carlos.mat q

% then run Model_Room.slx
