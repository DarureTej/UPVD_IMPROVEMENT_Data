clc
clear, clc, close all

Nbjour = 7; %Nombre de jours de simulation

Te=10; %P�riode d'�chantillonnage en minutes
t=(0:Te:Nbjour*24*Te*6-Te)'; %Vecteur temps

%Tarification sur une journ�e de 00h00 � 23h50
TEDF_jour(1:6*6,1) = 0.1228;
TEDF_jour(6*6+1:8*6,1) = 0.1579;
TEDF_jour(8*6+1:12*6,1) = 0.1228;
TEDF_jour(12*6+1:14*6,1) = 0.1579;
TEDF_jour(14*6+1:16*6,1) = 0.1228;
TEDF_jour(16*6+1:22*6,1) = 0.1579;
TEDF_jour(22*6+1:24*6,1) = 0.1228;

%Trac� du tarif de l'�lectricit� sur un jour
figure
plot(t(1:144)/60,TEDF_jour)
xlabel('Temps (h)')
ylabel('Prix (�/kWh)')
xlim([0 24])

%Tarification sur 7 jours
TEDF=[];
for i =1:Nbjour
    TEDF = [TEDF; TEDF_jour];
end

load data_exemple %Importation des donn�es de GHI et de consommation

%Calcul de la production solaire PV
S=40; %surface des panneaux
r = 0.20; %rendement des panneaux
Cp=0.1; %coefficient de pertes
Pp=S*r*GHI_real*(1-Cp); %production r�elle
Pp_pred=S*r*GHI_cs*(1-Cp); %production pr�dite

%Trac� pour comparaison des pr�diction et des valeurs r�elles
figure
subplot(3,1,1)
plot(temps_data,[GHI_real,GHI_cs])
ylabel('GHI (kW/m�)')
xlabel('Temps')
legend('GHI','GHI_{pred}')
subplot(3,1,2)
plot(temps_data,[Pp,Pp_pred])
ylabel('Pp (kW)')
xlabel('Temps')
legend('Pp','Pp_{pred}')
subplot(3,1,3)
plot(temps_data,[Pc Pc_pred])
ylabel('Puissance moyenne (kW)')
xlabel('Temps')
legend('Pc','Pc_{pred}')

FIT_Pc = 100*(1-norm(Pc-Pc_pred)/norm(Pc-mean(Pc)));
FIT_GHI = 100*(1-norm(GHI_real-GHI_cs)/norm(GHI_real-mean(GHI_real)));
FIT_Pp = 100*(1-norm(Pp-Pp_pred)/norm(Pp-mean(Pp)));

disp(['Fit de la pr�diction de la consommation : ',num2str(FIT_Pc),'%'])
disp(['Fit de la pr�diction du GHI : ',num2str(FIT_GHI),'%'])
disp(['Fit de la pr�diction de la production : ',num2str(FIT_Pp),'%'])


Capacite_batterie = 100;%Capacit� de la batterie en kWh
Eb(1) = 0.1*Capacite_batterie; %Initialisation de la batterie � 10% (valeur minimale)

%boucle de simulation et optimisation
Tfin = length(t); %Temps de simulation totale
Hp = 6*24; %Horizon de pr�diction en pas de temps
%Pr_opt=zeros(Hp,1);
for i=1:Tfin-1

    disp(' ')
    disp(['it�ration ',num2str(i),'/',num2str(Tfin-1)])
    %Optimisation de la commande
    
    Hpr = min([Tfin-i Hp]);% Horizon de pr�diction raccourci automatiquement � la fin de la simulation
    
    TEDF_p = TEDF(i+1:i+Hpr,1);%Pr�diction de la tarification 

    %Choix de la pr�diction utilis�e
    prediction = 0;
    switch prediction
        case 0
            Pc_p = Pc(i+1:i+Hpr,1);%Pr�diction parfaite de la consommation
            Pp_p = Pp(i+1:i+Hpr,1);%Pr�diction parfaite de la production
        case 1
            Pc_p = Pc(i+1:i+Hpr,1);%Pr�diction parfaite de la consommation
            Pp_p = Pp_pred(i+1:i+Hpr,1);%Pr�diction r�elle de la production
        case 2
            Pc_p = Pc_pred(i+1:i+Hpr,1);%Pr�diction r�elle de la consommation
            Pp_p = Pp(i+1:i+Hpr,1);%Pr�diction parfaite de la production
        case 3
            Pc_p = Pc_pred(i+1:i+Hpr,1);%Pr�diction r�elle de la consommation
            Pp_p = Pp_pred(i+1:i+Hpr,1);%Pr�diction r�elle de la production
    end
    
    %Diff�rentes conditions initiales pour l'optimisation
    Pr_0=ones(Hpr,1);
    %Pr_0=Pc_p-Pp_p;
    %Pr_0=[Pr_opt(2:Hpr); 1];
    
    Pr_min = -9*ones(size(Pr_0)); % - puissance maximale vendue kW
    Pr_max = 9*ones(size(Pr_0)); % puissance maximale achet�e kW
    
    %Choix de la m�thode d'optimisation
    methode = 4;
        
    switch methode
        case 0
            Pr_opt = Pr_0;
            fob = 0;
            exitflag = 0;
        case 1
            disp('methode = SQP')
            options = optimoptions(@fmincon,'Display','none','Algorithm','sqp','MaxFunEvals',10000,'MaxIter',1000);
            [Pr_opt,fob,exitflag] = fmincon(@(Pr_opt) funobj(Pr_opt,Eb(i,1),TEDF_p,Pc_p,Pp_p,Capacite_batterie,Te),Pr_0,[],[],[],[],Pr_min,Pr_max,@(Pr_opt) mycon(Pr_opt,Eb(i,1),Pc_p,Pp_p,Capacite_batterie,Pr_max,Pr_min,Te),options);
        case 2
            disp('methode = Recherche directe avec motifs g�n�ralis�s')
            options = psoptimset(@patternsearch);
            options = psoptimset(options,'Display','none');
            [Pr_opt,fob,exitflag] = patternsearch(@(Pr_opt) funobj(Pr_opt,Eb(i,1),TEDF_p,Pc_p,Pp_p,Capacite_batterie,Te),Pr_0,[],[],[],[],Pr_min,Pr_max,@(Pr_opt) mycon(Pr_opt,Eb(i,1),Pc_p,Pp_p,Capacite_batterie,Pr_max,Pr_min,Te),options);
        case 3
            disp('methode = Algorithme g�n�tique')
            options = gaoptimset(@ga);
            options = gaoptimset(options,'Display','none','MutationFcn',@mutationadaptfeasible,'Generations',1000*length(Pr_0));
            [Pr_opt,fob,exitflag] = ga(@(Pr_opt) funobj(Pr_opt,Eb(i,1),TEDF_p,Pc_p,Pp_p,Capacite_batterie,Te),length(Pr_0),[],[],[],[],Pr_min,Pr_max,@(Pr_opt) mycon(Pr_opt,Eb(i,1),Pc_p,Pp_p,Capacite_batterie,Pr_max,Pr_min,Te),options);
        case 4
            disp('methode = Essaims particulaires')
            fminconopts= optimoptions(@fmincon,'Display','none','Algorithm','sqp');
            options = optimoptions(@particleswarm,'Display','none','HybridFcn', {@fmincon,fminconopts});
            [Pr_opt,fob,exitflag] = particleswarm(@(Pr_opt) funobj2(Pr_opt,Eb(i,1),TEDF_p,Pc_p,Pp_p,Capacite_batterie,Pr_max,Pr_min,Te),length(Pr_0),Pr_min,Pr_max,options);
    end
    
    %Application du premier incr�ment de commande au syst�me
    Pr(i,1) = Pr_opt(1);
        
    Pb(i,1)=Pp(i,1)-Pc(i,1)+Pr(i,1); %Calcul de la puissance de stockage/d�stockage
    [Eb_s,Pb_real_s]=systeme(Pb(i,1),Eb(i,1),Capacite_batterie,Te); %Evolution de l'�tat de la batterie et puissance de stockage/d�stockage r�elle
    Pb_real(i,1) = Pb_real_s;
        
    Pr_real(i,1)=-Pp(i,1)+Pc(i,1)+Pb_real(i,1); %Calcul de la puissance r�ellement vendue/soutir�e
    
    if Pr_real(i,1)<Pr_min(1)||Pr_real(i,1)>Pr_max(1)
        disp('contraintes d�pass�es en r�el -----------------------------------------------')
    end
    
    Eb(i+1,1) = Eb_s;
    
    disp(['fob = ',num2str(fob)])
    disp(['exitflag = ',num2str(exitflag)])

end

Cr = (TEDF(1:end-1).*Pr_real*Te/60);%Calcul de l'�volution du co�t
Crtot = sum(TEDF(1:end-1).*Pr_real*Te/60);%Calcul du co�t total

Ebmin = Capacite_batterie*0.1*ones(size(Eb));
Ebmax = Capacite_batterie*ones(size(Eb));

Pr_min = -9*ones(size(Pr));
Pr_max = 9*ones(size(Pr));

%Trac� des graphiques
disp(' ')
disp('-------------------------------------------')
disp(['Co�t total : ',num2str(Crtot),' �']);

figure
subplot(3,1,1)
plot(t/60,Eb,t/60,Ebmin,':',t/60,Ebmax,':')
xlabel('Temps')
ylabel('Energie (kWh)')
legend('Eb','Ebmin','Ebmax')

subplot(3,1,2)
plot(t(1:end-1)/60,Pr_real,t(1:end-1)/60,Pb_real,t/60,Pc,'--',t/60,Pp,'--',t(1:end-1)/60,Pr_min,'b:',t(1:end-1)/60,Pr_max,'b:')
ylabel('Puissance kW')
xlabel('Temps')
legend('Pr','Pb','Pc','Pp','Prmin','Prmax')

subplot(3,1,3)
yyaxis left
plot(t(1:end-1)/60,Cr)
ylabel('Co�t (�)')
xlabel('Temps')
yyaxis right
plot(t/60,TEDF)
ylabel('Tarif (�/kWh)')
legend('Cr','Tarif')