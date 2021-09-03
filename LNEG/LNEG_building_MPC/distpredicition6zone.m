  function [distp,distt,distds]=distpredicition6zone(t,N)
load('dist')% loading disturbace
l=length(dist(:,1));
dist(:,1)=dist(:,1)+80*rand(l,1);
dist11=[dist(:,1)/1000/0.5 (dist(:,2) )]; % units in kW and °C 
dist1=dist11-[ones(l,1)*0.65 ones(l,1)*5];
dist1=[repmat(dist1(:,1),1,6) dist1(:,2)];
dist1=[dist1(:,1:4) dist1(:,5:6) dist1(:,7)];
q=dist1(t,1:6);Toa=dist1(t,7);
distt=[q';Toa];
d=distt;
for j=2:N
    d=[d  ,[dist1((t+j),1:6)';dist1((t+j),7)]];
end
distp=d;% for cetralized 
distds1=[dist1(:,1:3),dist1(:,7)];
dds=[dist1(t,1:3)';Toa]; % for distributed 
for j=2:N
    dds=[dds  ,[distds1((t+j),1:3)';distds1((t+j),4)]];
end
distds=dds;
