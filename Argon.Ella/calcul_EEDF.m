function [ne,eedf] = calcul_EEDF(ne_Te_filename,E)
%**********************************************************
% REceived from AMK JAnuary 2007
% Modified by Ella Sciamma (EMS) October 2007
% get the Te and ne flat profiles from file "ne_Te_filename.dat"
% calculates a Maxwellian electron energy distribution function
% from the flat Te and ne

%clear;

% number of energy points
% Esize=740;

% Initialize matrices to 0
eedf=zeros(length(E),1);
te=zeros;
ne=zeros;
% E=zeros(1,Esize);

% Define electron mass
me=9.11E-31;

%extra=1;

% Create E values
%E(1:11)=0:0.001:0.01;
%E(12:20)=0.02:0.01:0.1;
%E(21:29)=0.2:0.1:1;
%E(30:43)=1.2:0.2:3.8;
%E(44:Esize)=4:1:700;

extra=1;

% Get the flat Te and ne from file 'ne_Te_filename.dat'
% ne is measured experimentally, Te is obtained from Ar II ADAS CR model
ne_Te_dat=load(ne_Te_filename);
ne=ne_Te_dat(1);
Te=ne_Te_dat(2);

extra=1;

% Create a Maxwellian distribution from ne and Te
%for j=1:Esize
for j=1:length(E)
    if j==1 
        eedf(j)=ne*sqrt(me/(2*pi*1.6E-19*Te))*exp(-E(j)/Te)*(100*1.6E-19*(.001)/me);

    else
        eedf(j)=ne*sqrt(me/(2*pi*1.6E-19*Te))*exp(-E(j)/Te)*(100*1.6E-19*(E(j)-E(j-1))/me);%f(v)vdv to get cm^-2 s^-1
    end
end

% write the eedf's to a file that can be read into Bogaert's code
%save eedf_norad.dat eedf -ascii;

% write ne to a file that can be read into Bogaert's code
%save ne_norad.dat ne -ascii;
