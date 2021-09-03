function [fob] = funobj(Pr,Eb0,Tarif,Pc,Pp,Capacite_batterie,Te)
%Fonction de calcul de la fonction objectif

if isrow(Pr)
    Pr=Pr';
end

Pb=Pp-Pc+Pr; %Calcul de la puissance de stockage/déstockage

[~,Pb_real]=systeme(Pb,Eb0,Capacite_batterie,Te);%Calcul de la nouvelle capacité de la batterie et de la puissance réelle stockée/destockée

Pr_real=-Pp+Pc+Pb_real; %Calcul de la puissance réellement soutirée/vendue

fob = sum(Tarif.*Pr_real*Te/60);
end

