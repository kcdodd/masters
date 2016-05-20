function [eedf] = calcul_EEDF_3(Te ,E)
%**********************************************************
% REceived from AMK JAnuary 2007
% Modified by Ella Sciamma (EMS) October 2007
% Modified by K. Carter Dodd Nov. 2011, takes Te as argument, changed distribution to 3D isotropic


% Initialize matrices to 0
eedf=zeros(length(E),1);

% Define electron mass
me=9.11E-31;


% Create a Maxwellian distribution from ne and Te
%for j=1:Esize
for j=1:length(E)
    if j==1 
        %eedf(j)=sqrt(me/(2*pi*1.6E-19*Te))*exp(-E(j)/Te)*(100*1.6E-19*(.001)/me);
        eedf(j)=sqrt((2/pi)*(me/(1.6E-19*Te))^3)*(2*1.6E-19*E(j)/(me*me))*exp(-E(j)/Te)*(100*1.6E-19*(.001));
    else
        %eedf(j)=sqrt(me/(2*pi*1.6E-19*Te))*exp(-E(j)/Te)*(100*1.6E-19*(E(j)-E(j-1))/me);%f(v)vdv to get cm s^-1
        eedf(j)=sqrt((2/pi)*(me/(1.6E-19*Te))^3)*(2*1.6E-19*E(j)/(me*me))*exp(-E(j)/Te)*(100*1.6E-19*(E(j)-E(j-1)));
    end
end

