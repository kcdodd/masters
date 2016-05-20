function [rates] = arI_rates(Te_list)
% Te_list(eV): list of electron temperatures to compute solutions.
% ne_list (cm^-3): list of electron densities to compute solutions.
% n0_list (cm^-3): list of neutral densities to compute solutions.
% ion_fraction: fraction of electron density to use as ion density. (n_ion
% = ion_fraction*n_e)
% plasma_radius (cm): describes how big the plasma is to determin escape factors and meta-stable diffusion.
% allCombos: true = compute solution for every combination of Te, ne, and
% n0. false = only compute each set of {te, ne, n0}.
% start_level: first level to keep from solutions.
% end_level: last level to keep from solutions.
% return solutions: struct containing table of all precomputed level densities, and solution parameters.
% return varargout: details of each solution
%==========================================================
% Ar I Collisional radiative model (1D)
% received by Ella Sciamma(EMS)from Amy Keesee
% modified October 2007 by EMS 
% --> convertion from fortran to matlab
% --> took off radial dependency
% 
% get ne and Te from "ne_Te_filename" (now specified as aruments KCD)
% "factor" is used to calculate the neutral density from ne (neutral
% density is now explicitly specified KCD)
% --> if factor=0.25, n0=ne*0.25 => degree of ionization
%       ne/(ne+n0)=80%
% "xrad" is the radius of the plasma column (now plasma_radius KCD)
% --> used in the calculation of the metastable diffusion
% "press" is the fill pressure
% --> 0.05 mTorr for VX-100
% --> 1 mTorr for Helicon old and new and Helimak
%
%----------------------------------------------
% Modified by K. Carter Dodd (KCD), Oct. 2011
%----------------------------------------------
% Now computes a lookup table of solutions.
%
% Modified sequence of operations to be more efficient in
% computing multiple solutions. Only computes <sigma*v> once
% per designated electron temperature for all rates, and multplies by ne
% for each designated electron density. Finally, computes solutions
% for each designated neutral density.
%
% Added ability to change how many levels,
% and which processes, to include in the solutions.
%
% Changed function arguments and names.
%
% Computation of electron energy distribution function now in this file,
% Changed distribution to 3D isotropic, and now uses adaptive quadrature for integration
% instead of rectangle method. Had to vectorize cross-section functions.
%
% Corrected issue that caused error in calculation of escape factor in the limit
% of being optically thin.
%
% Corrected issue that interchaged excitation and de-excitation cross-sections between some levels.
%
% Removed variables that are not implemented.
%
% Changed iteration method and methods within the inner solution loops.
%
% Initial guess for solution is now LTE for iterative method instead of all in ground state.
%
% 
%==========================================================

useQuadGK = true; % set to false if quadgk is not supported

% how many of each parameter are there
numTe = size(Te_list, 2);

% INITIALIZATION - PARAMETERS
%============================
% Total number of level parents considered
Ntot=65;

% rate coefficient integrals without electron/ion density factors.
ion = zeros(numTe,Ntot);
recomb_2b = zeros(numTe,Ntot);
recomb_3b = zeros(numTe,Ntot);
trans_e = zeros(numTe,Ntot,Ntot);

% statistical weight factor
fac(Ntot)=zeros;

% Input file variables
Eion(Ntot)=zeros;
Eexc(Ntot)=zeros;
aAf(Ntot,Ntot)=zeros;
aP(Ntot,Ntot)=zeros;
gam(Ntot)=zeros;
g(Ntot)=zeros;
g0(Ntot)=zeros;
gi(Ntot)=zeros;
A(Ntot,Ntot)=zeros;


%parameter(pi=3.1415926536)
pi=3.1415926536;



%3. Input of data from files
%---------------------------
    % LEVELS1M.DAT: data necessary for cross sections
    %- - - - - - - - - - - - - - - - - - - - - - - - -
% read ionization+excitation energy, and degeneracies of the levels:
% the data are:
% n, Eion, Eexc, g(n), g0, gi, gam(n)
levels1m = load('levels1m.dat');
Eion=levels1m(:,2)'; %ATTENTION, Eion is a row, so we need to transpose the 2nd column of levels1m.dat
Eexc=levels1m(:,3)'; %same thing with Eexc
g=levels1m(:,4)'; % same thing with g
g0=levels1m(:,5)'; % same thing for g0
gi=levels1m(:,6)'; % same thing for gi
gam=levels1m(:,7)'; %same thing for gam

% Define fac(n),
% gi(n) is equal to 2, 4 or 6 and no other value
for n=1:Ntot
    if gi(n) == 6
        fac(n)=1.0;
    end
    if gi(n) == 4
        fac(n)=0.667;
    end
    if gi(n) == 2
        fac(n)=0.333;
    end
end

    %LEVELS2M.DAT: aA*fmn (allowed) and aP (forbidden) coeff.for elec.excit:
    %LEVELS3M.DAT: Amn (trans.probab.):
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% data in levels2m.dat are: n, m, aAf*f(n,m), aP(n,m)
% data in levels3m.dat are: n, m, Anm
levels2m = load('levels2m.dat');
levels3m = load('levels3m.dat');

for n=1:Ntot-1
    if n == 1
        low=n;
        high=Ntot-n;
    else
        low=low + Ntot - (n-1);
        high=high + Ntot - n;
    end
    aAf(n,n+1:Ntot)=levels2m(low:high,3)';
    aP(n,n+1:Ntot)=levels2m(low:high,4)';
    A(n,n+1:Ntot)=levels3m(low:high,3)';
end

trans_rad = A'; % want 1st index to mark source level.

for Te_i = 1:numTe

	disp(['Te= ' num2str(Te_list(Te_i))]); % track solutions progress
    

	%====================================================================
	% CALCULATION OF COLLISION 'RATES'
	%=================================


	% ELECTRONS
	%----------
	for n = 1:Ntot

		me = 9.11E-31; %electron mass is kg
		Te = Te_list(Te_i); % get electron temperature from list
        
        if (Te == 0)
            
            continue;
        end

		%energy ranges to use in integrations (eV). Integration is
		% broken up into ranges with limits E_i to E_j, starting at 
		% either E=0, or lowest energy.
        E_0 = 0.001;
		E_1 = 0.01;
		E_2 = 0.1;
		E_3 = 1;
		E_4 = 4;
		E_5 = 20;
		E_6 = 100;
		E_7 = 1000;

		% define the v*f part of sigma*v*f to compute
		% <sigma*v>: f=(3D isotropic maxwell-boltzmann). with v in cm/s, Te and Ee in eV.
		vf = @(Ee) 100*1.6E-19*sqrt((2/pi)*(me/(1.6E-19*Te))^3).*(2*1.6E-19*Ee/(me*me)).*exp(-Ee/Te);
		%vf = @(Ee) 100*1.6E-19*sqrt(me/(2*pi*1.6E-19*Te))*exp(-Ee/Te)/me; % 1-D version
		
		
		% define electron impact ionization integrand
		kioneU = @(Ee) vf(Ee).*Sione(Ee ,n,Eion(n)); % define sgima*v*f product
		
		% use adaptive quadrature for different energy ranges to capture
		% structure of cross-sections and distribution function.
		% integrate only E > Eion
		
		if useQuadGK
			ion(Te_i, n) = quadgk(kioneU, Eion(n), Eion(n)+E_1)+quadgk(kioneU, Eion(n)+E_1, Eion(n)+E_2)+quadgk(kioneU, Eion(n)+E_2, Eion(n)+E_3)+quadgk(kioneU, Eion(n)+E_3, Eion(n)+E_4)+quadgk(kioneU, Eion(n)+E_4, Eion(n)+E_5)+quadgk(kioneU, Eion(n)+E_5, Eion(n)+E_6)+quadgk(kioneU, Eion(n)+E_6, Eion(n)+E_7);
		else
			ion(Te_i, n) = quad(kioneU, Eion(n), Eion(n)+E_1)+quad(kioneU, Eion(n)+E_1, Eion(n)+E_2)+quad(kioneU, Eion(n)+E_2, Eion(n)+E_3)+quad(kioneU, Eion(n)+E_3, Eion(n)+E_4)+quad(kioneU, Eion(n)+E_4, Eion(n)+E_5)+quad(kioneU, Eion(n)+E_5, Eion(n)+E_6)+quad(kioneU, Eion(n)+E_6, Eion(n)+E_7);
		end
		
		
		% define radiative recombination integrand
		rrec1eU = @(Ee) vf(Ee).*Srec1e(Ee ,n,Eion(n),gam(n),g(n),gi(n)) * fac(n);
		
		if useQuadGK
			recomb_2b(Te_i, n) = quadgk(rrec1eU, E_0, E_1)+quadgk(rrec1eU, E_1, E_2)+quadgk(rrec1eU, E_2, E_3)+quadgk(rrec1eU, E_3, E_4)+quadgk(rrec1eU, E_4, E_5)+quadgk(rrec1eU, E_5, E_6)+quadgk(rrec1eU, E_6, E_7);
		else
			recomb_2b(Te_i, n) = quad(rrec1eU, 0, E_1)+quad(rrec1eU, E_1, E_2)+quad(rrec1eU, E_2, E_3)+quad(rrec1eU, E_3, E_4)+quad(rrec1eU, E_4, E_5)+quad(rrec1eU, E_5, E_6)+quad(rrec1eU, E_6, E_7);
		end
		
		% define three-body recombination integrand
		rrec2eU = @(Ee) vf(Ee).*Srec2e(Ee ,n,Eion(n),g(n),gi(n), Te) * fac(n);
		
		if useQuadGK
			recomb_3b(Te_i, n) = quadgk(rrec2eU, 0, E_1)+quadgk(rrec2eU, E_1, E_2)+quadgk(rrec2eU, E_2, E_3)+quadgk(rrec2eU, E_3, E_4)+quadgk(rrec2eU, E_4, E_5)+quadgk(rrec2eU, E_5, E_6)+quadgk(rrec2eU, E_6, E_7);
		else
			recomb_3b(Te_i, n) = quad(rrec2eU, 0, E_1)+quad(rrec2eU, E_1, E_2)+quad(rrec2eU, E_2, E_3)+quad(rrec2eU, E_3, E_4)+quad(rrec2eU, E_4, E_5)+quad(rrec2eU, E_5, E_6)+quad(rrec2eU, E_6, E_7);		    
		end

		for m = 1:Ntot
			% only transition to other levels
			if m ~= n
				% aAf and aP are only defined for k > j, so get correct values
				j = min(n,m);
				k = max(n,m);
			
				A_coeff = aAf(j, k);
				P_coeff = aP(j, k);
                    		
				% m > n does not mean Eexc(m) > Eexc(n), so determine if this is excitation or 
				% de-excitation
				if (Eexc(m) > Eexc(n))

					% this is an excitation transition
					Emn=Eexc(m)-Eexc(n);

					% define electron impact excitation integrand
					kexceU = @(Ee) vf(Ee).*Sexce(Ee ,n,m,Eexc(n),Eexc(m),A_coeff,P_coeff,g(n),g(m));

					%integrate only from E > Emn

					if useQuadGK
						trans_e(Te_i, n,m) = quadgk(kexceU, Emn, Emn+E_1)+quadgk(kexceU, Emn+E_1, Emn+E_2)+quadgk(kexceU, Emn+E_2, Emn+E_3)+quadgk(kexceU, Emn+E_3, Emn+E_4)+quadgk(kexceU, Emn+E_4, Emn+E_5)+quadgk(kexceU, Emn+E_5, Emn+E_6)+quadgk(kexceU, Emn+E_6, Emn+E_7);
					else
						trans_e(Te_i, n,m) = quad(kexceU, Emn, Emn+E_1)+quad(kexceU, Emn+E_1, Emn+E_2)+quad(kexceU, Emn+E_2, Emn+E_3)+quad(kexceU, Emn+E_3, Emn+E_4)+quad(kexceU, Emn+E_4, Emn+E_5)+quad(kexceU, Emn+E_5, Emn+E_6)+quad(kexceU, Emn+E_6, Emn+E_7);
					end

				else
					% this is a de-excitation transion

					% define electron impact de-excitation integrand
					kdeexeU = @(Ee) vf(Ee).*Sdeexe(Ee,n,m,Eexc(n),Eexc(m),A_coeff,P_coeff,g(n),g(m),g0(n));

					if useQuadGK
						trans_e(Te_i, n,m) = quadgk(kdeexeU, E_0, E_1)+quadgk(kdeexeU, E_1, E_2)+quadgk(kdeexeU, E_2, E_3)+quadgk(kdeexeU, E_3, E_4)+quadgk(kdeexeU, E_4, E_5)+quadgk(kdeexeU, E_5, E_6)+quadgk(kdeexeU, E_6, E_7);
					else
						trans_e(Te_i, n,m) = quad(kdeexeU, 0, E_1)+quad(kdeexeU, E_1, E_2)+quad(kdeexeU, E_2, E_3)+quad(kdeexeU, E_3, E_4)+quad(kdeexeU, E_4, E_5)+quad(kdeexeU, E_5, E_6)+quad(kdeexeU, E_6, E_7);
					end
				end

			end
            
		end
        
	end% n

end %Te_i



rates = struct('Te', Te_list, 'trans_e', trans_e, 'trans_rad', trans_rad, 'ion', ion, 'recomb_2b', recomb_2b, 'recomb_3b', recomb_3b);