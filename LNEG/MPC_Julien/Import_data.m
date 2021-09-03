clear, clc, close all

%----------------------------------------------------------------------
% Importation des données de consommation
%----------------------------------------------------------------------
%Consommation réelle
[num,txt,raw]=xlsread('Consommation.xlsx',1);

for i=0:6
    txt{2+i*48,1}=[txt{2+i*48,1},' 00:00:00'];
end
temps1 = datetime(txt(2:end,1),'InputFormat','dd/MM/yyyy HH:mm:SS');

pmoy1 = num(1:end-1,4)/1000;

t1 = temps1(1);
t2 = temps1(end);

temps2 = (t1:minutes(10):t2+minutes(20))';
pmoy2=interp1(temps1,pmoy1,temps2,'previous','extrap');


%Prédiction de consommation

num=xlsread('Consommation.xlsx',2,'H2:L337');
pmoy_moy1 = mean(num(:,3:5),2)/1000;
pmoy_moy2=interp1(temps1,pmoy_moy1,temps2,'previous','extrap');

FIT_Pc = 100*(1-norm(pmoy2-pmoy_moy2)/norm(pmoy2-mean(pmoy2)))

%Tracé prédiction de consommation

% figure
% stairs(temps1,[pmoy1 pmoy_moy1])
% ylabel('Puissance moyenne (kW)')
% xlabel('Temps')

figure
stairs(temps2,[pmoy2 pmoy_moy2])
ylabel('Puissance moyenne (kW)')
xlabel('Temps')
legend('Pc','Pc_{pred}')
% Importation des données de consommation

disp(['Conso1 = ',num2str(sum(pmoy1*30/60)),' kWh'])
disp(['Conso2 = ',num2str(sum(pmoy2*10/60)),' kWh'])

%----------------------------------------------------------------------
% Importation des données de GHI
%----------------------------------------------------------------------

num=xlsread('SoDa_HC3-METEO_lat42.681_lon2.891_2004-10-28_2004-11-04_1718803087.xlsx',1,'C33:D1040');
GHI_real = num(:,1)/(10/60)/1000;
GHI_cs = num(:,2)/(10/60)/1000;

FIT_GHI = 100*(1-norm(GHI_real-GHI_cs)/norm(GHI_real-mean(GHI_real)))

figure
stairs(temps2,[GHI_real,GHI_cs])
ylabel('GHI (kW/m²)')
xlabel('Temps')

temps_data = temps2;
Pc = pmoy2;
Pc_pred = pmoy_moy2;
save data_exemple temps_data Pc Pc_pred GHI_real GHI_cs
