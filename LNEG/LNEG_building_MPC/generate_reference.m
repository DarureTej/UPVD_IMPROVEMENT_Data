%% generate a refernce trajectory 
clc
load('dist')% loading disturbace
l=length(dist(:,1));
for i=1:l
if (dist(i,1)>0)
    Ulb(i)=-0.5;Uub(i)=0.5;
else
    Ulb(i)=-8.5;Uub(i)=8.5;
end
end 
% offset=ones(1,20);
% Ulb=[Ulb(1,5:end) Ulb(1,1:4)];Uub=[Uub(1,5:end) Uub(1,1:4)];
%  Ulb=circshift(Ulb,[1 6]);Uub=circshift(Uub,[1 6]);
% clearvars -except Ulb Uplot(Ulb)
