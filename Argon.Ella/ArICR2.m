function ArICR_table = ArICR2(ne_Te_filename,Nediv,experiment)
%====================================================================
% by Ella Sciamma
% October 2007
% 1. Define the experiment's plasma radius "xrad" and pressure "pres"
% 2. use Ar I CR model to calculate the Ar I populations
% of lvels 7, 8, 9, 10, and 11 for the neutral densities
% corresponding to the degrees of ionization:
% 1% 10% 20% 30% 40% 50% 60% 70% 80% 90% 99%
% save the Ar I populations of interest for all n0 in a unique file
%====================================================================
switch experiment
    case 1 %VX-100
        xrad=8; % length of plasma observed after ICRH=16 cm
        %pres=7.5; %mTorr
    case 2 %Helimak
        xrad=50; % length of plasma observed = 200 cm
        % but shortest length of plasma to reach the walls from
        % the plasma center =50 cm=radial length/2
        %pres=0.018; %mTorr
    case 3 % Helicon old
        xrad=2; % length of plasma observed after ICRH=4 cm (off centered)
        %pres=1 %mTorr
    case 4 % Helicon new
        xrad=3; % length of plasma observed after ICRH=6 cm
        %pres=1 %mTorr
end

% ArI CR population calculation for
%  1% degree of ionization
% n0=99*ne
[n0_1percent, ArI_pop_1percent] = arcrmodel_new2(ne_Te_filename,99, Nediv, xrad);

% 10% degree of ionization
% n0=9*ne
[n0_10percent, ArI_pop_10percent] = arcrmodel_new2(ne_Te_filename,9, Nediv, xrad);

% 20% degree of ionization
% n0=4*ne
%[n0_20percent, ArI_pop_20percent] = arcrmodel_new2(ne_Te_filename,4, Nediv, xrad);

% 30% degree of ionization
% n0=7/3*ne;
%[n0_30percent, ArI_pop_30percent] = arcrmodel_new2(ne_Te_filename,2.3333, Nediv, xrad);

% 40% degree of ionization
% n0=3/2*ne
%[n0_40percent, ArI_pop_40percent] = arcrmodel_new2(ne_Te_filename,1.5, Nediv, xrad);

% 50% degree of ionization
% n0=ne
%[n0_50percent, ArI_pop_50percent] = arcrmodel_new2(ne_Te_filename,1, Nediv, xrad);

% 60% degree of ionization
% n0=2/3*ne
%[n0_60percent, ArI_pop_60percent] = arcrmodel_new2(ne_Te_filename,0.6667, Nediv, xrad);

% 70% degree of ionization
% n0=3/7*ne
%[n0_70percent, ArI_pop_70percent] = arcrmodel_new2(ne_Te_filename,0.4286, Nediv, xrad);

% 80% degree of ionization
% n0=1/4*ne
%[n0_80percent, ArI_pop_80percent] = arcrmodel_new2(ne_Te_filename,0.25, Nediv, xrad);

% 90% degree of ionization
% n0=1/9*ne
%[n0_90percent, ArI_pop_90percent] = arcrmodel_new2(ne_Te_filename,0.1111, Nediv, xrad);

% 99% degree of ionization
% n0=1/99*ne
%[n0_99percent, ArI_pop_99percent] = arcrmodel_new2(ne_Te_filename,0.0101, Nediv, xrad);

levels(1:65, 1)=1:1:65;
%allpop=[levels ArI_pop_1percent ArI_pop_10percent ArI_pop_20percent ArI_pop_30percent ArI_pop_40percent ArI_pop_50percent ArI_pop_60percent ArI_pop_70percent ArI_pop_80percent ArI_pop_90percent ArI_pop_99percent];
allpop=[levels ArI_pop_1percent ArI_pop_10percent];
%ArICR_table(1,:)=[0 1 10 20 30 40 50 60 70 80 90 99];
ArICR_table(1,:)=[0 1 10 ];
%ArICR_table(2,:)=[0 n0_1percent n0_10percent n0_20percent n0_30percent n0_40percent n0_50percent n0_60percent n0_70percent n0_80percent n0_90percent n0_99percent];
ArICR_table(2,:)=[0 n0_1percent n0_10percent];
ArICR_table(3:6,:)=allpop(8:11,:);
