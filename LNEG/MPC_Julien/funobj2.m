function [fob] = funobj2(Pr,Eb0,Tarif,Pc,Pp,Capacite_batterie,Pr_max,Pr_min,Te)
%Fonction de calcul de la fonction objectif

if isrow(Pr)
    Pr=Pr';
end

Pb=Pp-Pc+Pr; %Calcul de la puissance de stockage/d�stockage

[~,Pb_real]=systeme(Pb,Eb0,Capacite_batterie,Te);%Calcul de la nouvelle capacit� de la batterie et de la puissance r�elle stock�e/destock�e

Pr_real=-Pp+Pc+Pb_real; %Calcul de la puissance r�ellement soutir�e/vendue

%Contraintes molles
lambda = 1000000;
Cm=lambda*((Pr_real>Pr_max).*(Pr_real-Pr_max)+(Pr_real<Pr_min).*(Pr_min-Pr_real));

fob = sum(Tarif.*Pr_real*Te/60) + sum(Cm);
end

