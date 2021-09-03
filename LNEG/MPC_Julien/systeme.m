function [Eb,Pb_real] = systeme(Pb,Eb0,Capacite_batterie,Te)
%Fonction de simulation du modèle non linéaire
Ebmin = 0.1*Capacite_batterie;
Ebmax = 1*Capacite_batterie;
Eb(1)=Eb0;
Pb_real = Pb;
for i=1:length(Pb)
    Eb(i+1,1)=Eb(i)+Pb(i)*Te/60; %Charge/décharge de la batterie
    
    %saturation de la capacité de la batterie et calcul de la puissance
    %réelle de stockage / déstockage
    if Eb(i+1,1)>Ebmax
        Eb(i+1,1)=Ebmax;
        Pb_real(i) = (Ebmax-Eb(i,1))*60/Te;
    elseif Eb(i+1,1)<Ebmin
        Eb(i+1,1)=Ebmin;
        Pb_real(i) = (Ebmin-Eb(i,1))*60/Te;
    end
end
Eb=Eb(2:end);

end

