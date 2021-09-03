% Generation of A, B, C, D matrices for the Lisboa meteo station
% Model parameters - 5R1C network model of EN ISO 13790
% calculated according to EN ISO 13790
cm = 37582341/3600 % J/K
am = 437.44 % m2
atot = 661.94 % m2
af = 156 % m2

htrop = 80.64 % W/K
htrw = 19.05 % W/K
htris = 3.45*atot
%hve - varying, computed by this program
htrms = 9.1*am
htrem = 1/(1/htrop-1/htrms)

% == A, B, C and D matrices and U vector for the whole year - 8760 hours ==
% weather data from Lisboa.xls file, days of the week numbered from 1 to 7 %
% Solar irradiance from EnergyPlus - Perez model

timehour = xlsread ('tmyLisboa', 1, 'a2:a8761') % hourly time step in table
timehour = timehour' % single-line vector 

month = xlsread ('tmyLisboa', 1, 'b2:b8761') % consecutive months
month = month' % single-line vector 

day = xlsread ('tmyLisboa', 1, 'c2:c8761') % days of a week (1-7)
day = day' % single-line vector 

% vector of the consecutive hours
hour = xlsread ('tmyLisboa', 1, 'd2:d8761') % hours from xls file
hour = hour' % single-line vector 


% hours according to legal time
lth = xlsread ('tmyLisboa', 1, 'm2:m8761') % legal time hours
lth = lth' % single-line vector 


% consecutive hours of a year (number)
wymnumber = size(timehour)
wymnumber = wymnumber(1,2)

% Time-dependent ventilation and internal gains schedules

for i = 1:wymnumber

switch day(i)
case {1, 2, 3, 4, 5} % weekdays (Mon-Fri) 
        switch lth(i)
            case {6, 7, 15, 16, 17, 18, 19, 20, 21} % hours 6-8, 15-22 
            hve(i) = 122.67 % 1 air change per hour (1 ACH)
            qint(i) = 800 % convective heat gains 
            
            case {8, 9, 10, 11, 12, 13, 14} % hours 8-15 
            hve(i) = 0 % 0 ACH
            qint(i) = 0 % convective heat gains 
                                    
            otherwise % night hours: 22-6
            hve(i) = 61.33 % 0,5 ACH
            qint(i) = 300 % convective heat gains 
        end

        
             
        
    otherwise % Sat-Sun
       
hve(i) = 122.67 % 1 ACH
qint(i) = 400
end
qqint=qint(i)
end




%f = @(x) x+2
%u=1:7
%arrayfun(f,u)

% denominator = (htrms + htrw)*(htris + hve) + hve*htris

% Computation of the A matrix, parameter hve as the function's argument
% anonymous function
aa11 = @(hve) (((htris+hve)*htrms*htrms/((htrms + htrw)*(htris + hve) + hve*htris))-htrem-htrms)/cm
%x=hhve
%u=1:7
a11=arrayfun(aa11,hve) % a11 element of the A matrix



% Computation of the B matrix, parameter hve as the function's argument
% anonymous function
% b [1x6]
bb11 = @(hve) (htrem+(htrms*htrw*(htris+hve))/((htrms + htrw)*(htris + hve) + hve*htris))/cm %Te
b11=arrayfun(bb11,hve) % element b11 of the matrix B

bb12 = @(hve) (htrms*hve*htris)/((htrms + htrw)*(htris + hve) + hve*htris)/cm %Tsup
b12=arrayfun(bb12,hve) % element b12 

bb13 = @(hve) 1/cm  % qm
b13=arrayfun(bb13,hve) % element b13 

bb14 = @(hve) ((htris+hve)*htrms/((htrms + htrw)*(htris + hve) + hve*htris))/cm %qst
b14=arrayfun(bb14,hve) % czyli element b14

bb15 = @(hve) htrms*htris/((htrms + htrw)*(htris + hve) + hve*htris)/cm %qia
b15=arrayfun(bb15,hve) % element b15

bb16 = @(hve) htrms*htris/((htrms + htrw)*(htris + hve) + hve*htris)/cm %qHC
b16=arrayfun(bb16,hve) % element b16


% Computation of the C matrix, parameter hve as the function's argument
% anonymous function
% c [2x1] -> [Ti; Ts]
cc11 = @(hve) htrms*htris/((htrms + htrw)*(htris + hve) + hve*htris) %Ti - Tindoor (internal temperature)
c11=arrayfun(cc11,hve) % element c11

cc21 = @(hve) htrms*(hve+htris)/((htrms + htrw)*(htris + hve) + hve*htris) %Ts - Tsurface
c21=arrayfun(cc21,hve) % element c21



% Computation of the D matrix, parameter hve as the function's argument
% anonymous function
% d [2x6] -> [Ti; Ts]

dd11 = @(hve)  htris*htrw/((htrms + htrw)*(htris + hve) + hve*htris) %Te
d11=arrayfun(dd11,hve) % element d11

dd12 = @(hve) ((htrms+htris+htrw)*hve)/((htrms + htrw)*(htris + hve) + hve*htris) %Tsup
d12=arrayfun(dd12,hve) % element d12

dd13 = @(hve)  0 % qm
d13=arrayfun(dd13,hve) % element d13 

dd14 = @(hve)  htris/((htrms + htrw)*(htris + hve) + hve*htris) %qst
d14=arrayfun(dd14,hve) % element d14 

dd15 = @(hve)  (htrms+htris+htrw)/((htrms + htrw)*(htris + hve) + hve*htris) %qia
d15=arrayfun(dd15,hve) % element d15 

dd16 = @(hve)  (htrms+htris+htrw)/((htrms + htrw)*(htris + hve) + hve*htris) %qHC
d16=arrayfun(dd16,hve) % element d16 

dd21 = @(hve)  (hve+htris)*htrw/((htrms + htrw)*(htris + hve) + hve*htris) %Te
d21=arrayfun(dd21,hve) % element d21 

dd22 = @(hve)  (htris*hve)/((htrms + htrw)*(htris + hve) + hve*htris) %Tsup
d22=arrayfun(dd22,hve) % element d22 

dd23 = @(hve)  0 %qm
d23=arrayfun(dd23,hve) % element d23

dd24 = @(hve)  (hve + htris)/((htrms + htrw)*(htris + hve) + hve*htris) %qst 
d24=arrayfun(dd24,hve) % element d24

dd25 = @(hve)  htris/((htrms + htrw)*(htris + hve) + hve*htris) %qia
d25=arrayfun(dd25,hve) % element d25

dd26 = @(hve)  htris/((htrms + htrw)*(htris + hve) + hve*htris) %qHC
d26=arrayfun(dd26,hve) % element d26




% Matrix A
A=[timehour; a11]
save ('ALisboa.mat','A')


%B=[b11,b12,b13,b14,b15,b16]
B=[timehour;b11;b12;b13;b14;b15;b16]
save ('BLisboa.mat','B')


% C=[c11; c21]

C=[timehour; c11; c21]
save ('CLisboa.mat','C')

% Matrix D converted to the form readable for the stvmgain function
D=[timehour; d11; d12; d13; d14; d15; d16; d21; d22; d23; d24; d25; d26]
save ('DLisboa.mat','D')

% Matrix D in normal form
Dstv = [timehour; d11; d21; d12; d22; d13; d23; d14; d24; d15; d25; d16; d26]
save ('DstvLisboa.mat','Dstv')




% Effective solar area of external surfaces Asol - from EN ISO 13790
aen = 0.2331 % N
aes = 3.6218 % S
aew = 2.2613 % W
aee = 2.1377 % E

aedw = 0.665 % Roof W
aede = 0.665 % Roof E

% coefficient for long-wave irradiance to sky - vertical walls
rsky = 3.6416

% coefficient for long-wave irradiance to sky - roof
rdsky = 1.876256

% Air temperature and sky temperature
te = xlsread ('tmyLisboa', 1, 'e2:e8761') % DBT - Dry Bulb Temperature
tsky = xlsread ('tmyLisboa', 1, 'f2:f8761') % TSKY - Sky temperature


% Global solar irradiance from TMY
qsn90 = xlsread ('tmyLisboa', 1, 'i2:i8761') % North
qss90 = xlsread ('tmyLisboa', 1, 'k2:k8761') % South
qsw90 = xlsread ('tmyLisboa', 1, 'l2:l8761') % West
qse90 = xlsread ('tmyLisboa', 1, 'j2:j8761') % East
qsdw = xlsread ('tmyLisboa', 1, 'h2:h8761') % West45
qsde = xlsread ('tmyLisboa', 1, 'g2:g8761') % East45


% Solar gains:
qqsn = qsn90*aen
qqss = qss90*aes
qqsw = qsw90*aew
qqse = qse90*aee

% Solar gains - roof:
qqd = qsdw*aedw+ qsde*aede

btr = 0.02651 % 1-btr factor - roof

%Long-wave gains - walls:
qlw = rsky*(te-tsky)

%Long-wave gains - roof
qdlw = rdsky*(te-tsky)

% Total solar
qsol = qqsn + qqss + qqsw + qqse + qqd*btr - qlw - qdlw*btr
qsol = qsol' % Total solar gains - vector


%================================================================


% Total internal heat gains qint as input file for Matlab
% heat fluxes to the 5R1C model nodes:
qia = 0.5*qint
qm = (am/atot)*(0.5*qint + qsol)
qst = (1-(am/atot)-(htrw/(9.1*atot)))*(0.5*qint + qsol)



% save ('qint1_10','qqint') % reserve

te=te'
q=[timehour; te; te; qm; qst; qia] % heat sources matrix (inputs) , 
% twice te because tsupply = te -> see EN ISO 13790
% q = [hour; Te; Tsup; qm; qst; qia] + qHC

save ('uLisboa.mat','q')

