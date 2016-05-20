function [solutions] = arcrmodel_2011(Te_list, ne_list, n0_list, ion_fraction, plasma_radius, allCombos, start_level, end_level)
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

% Which things to include in model
Nuse=65; % how many levels to include: 1 up to max is 65.
useIonRecomb = 1; % use ionization and recombination. =1 to use process, 0 to not use it.
useRadDecay = 1; % use spontaneous radiative decay
useThermal = 1; % Use collisions with thermal atoms. 
useDiff = 0; % use diffusion of meta-stable levels
useMeta = 0; % use losses between meta-stable levels and ground.
useQuadGK = false; % set to false if quadgk is not supported
maxIterations = 100000; % maximum number of iterations before giving up on solution


% how many of each parameter are there
numTe = size(Te_list, 2);
numNe = size(ne_list, 2);
numN0 = size(n0_list, 2);

if (~allCombos && (numTe ~= numNe || numTe ~=numN0))
	% if not computing all combonations then must have same number of
	% electron temperature, density and neutral density parameters.
	disp ('lists of pairs must have same length');
	return;
end

if (end_level < start_level)
	disp('end level must be higher then start');
end

if (end_level > Nuse)
	disp('cannot compute levels not included in model');
end

% how many levels from each solutions will be kept.
num_levels = end_level - start_level + 1;

% initialize total solution matrix
if allCombos
	level_table = zeros(num_levels, numTe, numNe, numN0);
else
	level_table = zeros(num_levels, numTe);
end

% INITIALIZATION - PARAMETERS
%============================
% Total number of level parents considered
Ntot=65;

% Number of energy point considered for EEDF
NE=739;

% Population variable
pop(Nuse,1)=zeros;

% Metastable diffusion, 2- and 3-body recombination coefficient rates
D(4)=zeros;
k2b(4)=zeros;
k3b(4)=zeros;

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

% old populations from previous time step
popold(Nuse,1)=zeros;

% Escape factor
Esc(Ntot,1)=zeros;

% Electron energy for EEDF
Ee(NE+1)=zeros;



% metastable variables
bbi(4,1)=zeros;
beta=zeros;
gamma=zeros;


%parameter(pi=3.1415926536)
pi=3.1415926536;

% Level parent indexes
n2(Ntot)=zeros;

%=================================================================
% DATA INPUT
%============
%1. General: ambient gas temperature and time step
%-------------------------------------------------------
% xtgas is the ambient gas temperature in K 
xtgas=348;
pres = 1;
dt=1e-3; % time steps

%2. Data for metastables
%------------------------
    % Diffusion coefficient to account for loss of metastables
    % D=74.6 cm^2 s^-1 at 1 Torr and 300 K
    % => linear equation to calculate D using press and xtgas
D(2)=74.6 * 1000/pres * xtgas/300;
D(4)=74.6 * 1000/pres * xtgas/300;
%dis('D(2)=')
%disp(D(2))
%disp('D(4)=')
%disp(D(4))

    % metastable-metastable atom collision rate coeff
xkmet=6.4e-10;

    % 2- and 3-body recombination with thermal ground state atoms
k2b(2)=2.3e-15;
k3b(2)=1.4e-32;
k2b(4)=4.3e-15;
k3b(4)=1.5e-32;


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


for Te_i = 1:numTe

	disp(['Te= ' num2str(Te_list(Te_i))]); % track solutions progress
    

	%====================================================================
	% CALCULATION OF COLLISION 'RATES'
	%=================================

	% rate coefficient integrals without electron/ion density factors.
	kioneIntegral=zeros(Ntot,1);
	rrec1eIntegral=zeros(Ntot,1);
	rrec2eIntegral=zeros(Ntot,1);
	kexceIntegral=zeros(Ntot,Ntot,1);
	kdeexeIntegral=zeros(Ntot,Ntot,1);


	% ELECTRONS
	%----------
	for n = 1:Ntot

		me = 9.11E-31; %electron mass is kg
		Te = Te_list(Te_i); % get electron temperature from list

		%energy ranges to use in integrations (eV). Integration is
		% broken up into ranges with limits E_i to E_j, starting at 
		% either E=0, or lowest energy.
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
			kioneIntegral(n) = quadgk(kioneU, Eion(n), Eion(n)+E_1)+quadgk(kioneU, Eion(n)+E_1, Eion(n)+E_2)+quadgk(kioneU, Eion(n)+E_2, Eion(n)+E_3)+quadgk(kioneU, Eion(n)+E_3, Eion(n)+E_4)+quadgk(kioneU, Eion(n)+E_4, Eion(n)+E_5)+quadgk(kioneU, Eion(n)+E_5, Eion(n)+E_6)+quadgk(kioneU, Eion(n)+E_6, Eion(n)+E_7);
		else
			kioneIntegral(n) = quad(kioneU, Eion(n), Eion(n)+E_1)+quad(kioneU, Eion(n)+E_1, Eion(n)+E_2)+quad(kioneU, Eion(n)+E_2, Eion(n)+E_3)+quad(kioneU, Eion(n)+E_3, Eion(n)+E_4)+quad(kioneU, Eion(n)+E_4, Eion(n)+E_5)+quad(kioneU, Eion(n)+E_5, Eion(n)+E_6)+quad(kioneU, Eion(n)+E_6, Eion(n)+E_7);
		end
		
		
		% define radiative recombination integrand
		rrec1eU = @(Ee) vf(Ee).*Srec1e(Ee ,n,Eion(n),gam(n),g(n),gi(n)) * fac(n);
		
		if useQuadGK
			rrec1eIntegral(n) = quadgk(rrec1eU, 0, E_1)+quadgk(rrec1eU, E_1, E_2)+quadgk(rrec1eU, E_2, E_3)+quadgk(rrec1eU, E_3, E_4)+quadgk(rrec1eU, E_4, E_5)+quadgk(rrec1eU, E_5, E_6)+quadgk(rrec1eU, E_6, E_7);
		else
			rrec1eIntegral(n) = quad(rrec1eU, 0, E_1)+quad(rrec1eU, E_1, E_2)+quad(rrec1eU, E_2, E_3)+quad(rrec1eU, E_3, E_4)+quad(rrec1eU, E_4, E_5)+quad(rrec1eU, E_5, E_6)+quad(rrec1eU, E_6, E_7);
		end
		
		% define three-body recombination integrand
		rrec2eU = @(Ee) vf(Ee).*Srec2e(Ee ,n,Eion(n),g(n),gi(n), Te) * fac(n);
		
		if useQuadGK
			rrec2eIntegral(n) = quadgk(rrec2eU, 0, E_1)+quadgk(rrec2eU, E_1, E_2)+quadgk(rrec2eU, E_2, E_3)+quadgk(rrec2eU, E_3, E_4)+quadgk(rrec2eU, E_4, E_5)+quadgk(rrec2eU, E_5, E_6)+quadgk(rrec2eU, E_6, E_7);
		else
			rrec2eIntegral(n) = quad(rrec2eU, 0, E_1)+quad(rrec2eU, E_1, E_2)+quad(rrec2eU, E_2, E_3)+quad(rrec2eU, E_3, E_4)+quad(rrec2eU, E_4, E_5)+quad(rrec2eU, E_5, E_6)+quad(rrec2eU, E_6, E_7);		    
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
						kexceIntegral(n,m)= quadgk(kexceU, Emn, Emn+E_1)+quadgk(kexceU, Emn+E_1, Emn+E_2)+quadgk(kexceU, Emn+E_2, Emn+E_3)+quadgk(kexceU, Emn+E_3, Emn+E_4)+quadgk(kexceU, Emn+E_4, Emn+E_5)+quadgk(kexceU, Emn+E_5, Emn+E_6)+quadgk(kexceU, Emn+E_6, Emn+E_7);
					else
						kexceIntegral(n,m)= quad(kexceU, Emn, Emn+E_1)+quad(kexceU, Emn+E_1, Emn+E_2)+quad(kexceU, Emn+E_2, Emn+E_3)+quad(kexceU, Emn+E_3, Emn+E_4)+quad(kexceU, Emn+E_4, Emn+E_5)+quad(kexceU, Emn+E_5, Emn+E_6)+quad(kexceU, Emn+E_6, Emn+E_7);
					end

				else
					% this is a de-excitation transion

					% define electron impact de-excitation integrand
					kdeexeU = @(Ee) vf(Ee).*Sdeexe(Ee,n,m,Eexc(n),Eexc(m),A_coeff,P_coeff,g(n),g(m),g0(n));

					if useQuadGK
						kdeexeIntegral(n,m) = quadgk(kdeexeU, 0, E_1)+quadgk(kdeexeU, E_1, E_2)+quadgk(kdeexeU, E_2, E_3)+quadgk(kdeexeU, E_3, E_4)+quadgk(kdeexeU, E_4, E_5)+quadgk(kdeexeU, E_5, E_6)+quadgk(kdeexeU, E_6, E_7);
					else
						kdeexeIntegral(n,m) = quad(kdeexeU, 0, E_1)+quad(kdeexeU, E_1, E_2)+quad(kdeexeU, E_2, E_3)+quad(kdeexeU, E_3, E_4)+quad(kdeexeU, E_4, E_5)+quad(kdeexeU, E_5, E_6)+quad(kdeexeU, E_6, E_7);
					end
				end

			end
            
		end
        
	end% n
    

	% make list of electron density indecies
	if allCombos
        	% compute all combonations of electron density and temperature.
        	computeNes = 1:size(ne_list,2);
	else
        	% only compute pairs of electron density and temperature.
        	computeNes = [Te_i Te_i];
	end % allCombos
        
	for ne_i = computeNes

		% get current electron density
		nes = ne_list(ne_i);
        
        	%disp(['      ne= ' num2str(nes)])

		% Define ion density. It can be different from ne if Nediv~=1
		ni=nes*ion_fraction;

		% multiply all integrals by factors of ni and ne. now k**** has units of cm^3/s, and r**** units of 1/s.
		kione = nes*kioneIntegral;
		rrec1e = ni*nes*rrec1eIntegral;
		rrec2e = ni*nes*nes*rrec2eIntegral;
		kexce = nes*kexceIntegral;
		kdeexe = nes*kdeexeIntegral;
        
		% make list of neutral density indecies
		if allCombos
			% compute all combonations of electron density and temperature.
			computeN0s = 1:size(n0_list,2);
		else
			% only compute pairs of electron density, temperature, and neutral density.
			computeN0s = [Te_i Te_i];
		end % allCombos
        
		for n0_i = computeN0s

			% current neutral density
			xn0 = n0_list(n0_i);
            
            		%disp(['            n0= ' num2str(xn0)])


			kionth=zeros(Ntot,1);
			rrecth=zeros(Ntot,1);
			kexcth=zeros(Ntot,Ntot);
			kdeexth=zeros(Ntot,Ntot);
			Kloss=zeros(Ntot,1);
			Kgain=zeros(Ntot,Ntot);
			Rgain=zeros(Ntot,1);

			% THERMALIZED ATOMS (E=0.03 eV -> v=3.81e4 cm/s) (k: per pop=1, in s-1)
			%-------------------
			Eth=0.03; % thermal energy
			vth=3.81e4; % thermal velocity

			for n = 1:Ntot
				kionth(n)=Sionth(Eion(n)) * vth * xn0;

				rrecth(n)=Sreci(Eth,n,Eion(n),g(n),gi(n)) * vth * xn0 * nes * fac(n) * ni;

				for m = 1:Ntot
					if m ~= n
						if (Eexc(m) > Eexc(n))

							% excitation
							kexcth(n,m)=Sexcth(n,m,Eexc(n),Eexc(m),gi(n),gi(m)) * vth * xn0;
						else
							% de-excitation
							kdeexth(n,m)=Sdeexth(n,m,Eexc(n),Eexc(m),g(n),g(m),gi(n),gi(m),g0(n)) * vth * xn0;
						end

					end

				end

			end% n

			%==========================================================================
			% CALCULATION OF THE ESCAPE FACTORS
			%==================================
			for n = 2:Ntot
				if (A(1,n) ~= 0.0)
					xkR=(2.1e-17)* g(n) * A(1,n) * xn0 * plasma_radius/((Eexc(n)^3)*sqrt(xtgas));
					xa=A(1,n) * (1 + (3.225e-14)* g(n) * xn0/(Eexc(n)^3))*(4.839e-9)/(Eexc(n)*sqrt(xtgas));
					Tc=sqrt(xa/(sqrt(pi)*xkR));

					if (xkR > 1.0)
						% xkR does not cause imaginary numbers if > 1.
						Td=1/(xkR*sqrt(pi*log(xkR)));
						Tcd=2*xa / (pi*sqrt(log(xkR)));

						Esc(n) = 1.9*Td*exp(-pi*(Tcd^2)/(4*Tc^2))+1.3*Tc*erf(sqrt(pi)*Tcd/(2*Tc));
					else
						% use limit as xkR -> 1 so that escape factor has
						% continuous value from xkR > 1 to xkR < 1.
						Esc(n) = 1.3*Tc;
					end

					if (Esc(n) > 1.0)
						% if for some reason it becomes bigger than one, create
						% explicit limit.
						Esc(n) = 1.0;
					end
				else
					Esc(n)=1.0;
				end

			end % n

			%==================================================================
			% CALCULATION OF THE LEVEL POPULATIONS
			%======================================
            

			% pre-compute rate terms for all levels which do not depend on populations:
			%----------------------------------------------------
			for n = 2:Nuse

				% electron, ion, fast atom, therm.atom ionization
				%- - - - - - - - - - - - - - - - - - - - - - - - - 
				Kloss(n) = useIonRecomb*kione(n) + useThermal*kionth(n);

				if (n == 2 || n == 4)
					Kloss(n)=Kloss(n) + useDiff*(D(n)/plasma_radius * (2/plasma_radius));
				end

				% radiative, 3-body, and thermal recombination
				%- - - - - - - - - - - - - - - - - - - - - - - - - 
				Rgain(n) = useIonRecomb*(rrec1e(n) + rrec2e(n) + useThermal*rrecth(n));

				for m = 1:Nuse


					if m ~= n
						if (Eexc(m) > Eexc(n))

							% excitation

							% electron, ion, fast atom, therm.atom excitation to higher levels
							%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

							Kloss(n) = Kloss(n) + kexce(n,m) + useThermal*kexcth(n,m);

							Kgain(n,m) = kdeexe(m,n) + useThermal*kdeexth(m,n) + useRadDecay*A(n,m);
						else
							%de-excitation

							% elec,ion,fast atom,therm.atom deexcit, radiat.decay to lower levels
							%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
							% (escape factors to level 1 incorporated)

							if m == 1
								Arad=A(m,n) * Esc(n);
							end

							if m ~= 1
								Arad=A(m,n);
							end

							Kloss(n)= Kloss(n) + kdeexe(n,m) + useThermal*kdeexth(n,m) + useRadDecay*Arad;

							Kgain(n,m) = kexce(m,n) + useThermal*kexcth(m,n);
						end
					end

				end % m
			    
			    
			end % n




			% initial guess for solution is LTE
			pop(1) = 1;
	    
			for n = 2:Nuse
				pop(n) = (g(n)/g(1))*exp(-Eexc(n)/Te);
			end

			% normalize to current neutral density to statisfy condition sum(pop) == xn0.
			pop = (xn0/sum(pop))*pop;

			it = 0;

			
			%solve for the solution iterativly. if max number is reached no solution convergence
			for it = 1:maxIterations
			    
				% old values
				%- - - - - - -
				popold=pop; %vectorized loop

				% compute iteration for all levels
				%+++++++++++++
				for n = 2:Nuse

					meta_loss = 0.0;

					if (useMeta == 1 && (n == 2 || n == 4))

						% additional losses for meta-stable levels (meta-meta collisions + meta-ground collisions)
						%...............................................
						if n == 2
						    pop2=pop(4);
						end
	
						if n == 4
						    pop2=pop(2);
						end
	
						meta_loss = 2*xkmet*pop(n) + xkmet*pop2 + pop(1)*k2b(n) + pop(1)^2*k3b(n);
					end

					% production from other levels and recombination
					% Kgain(n,:)*pop = row vector times column vector (inner product)
					% as loop continues uses progession of solution to calculate rates.
					%..................................
					prod = Kgain(n,:)*pop + Rgain(n);

					% iterate:
					% solution defined as d(pop(n))/dt = prod - pop(n)*(meta_loss + Kloss(n)) = 0.
					% gives: pop(n) = prod/(meta_loss + Kloss(n))
					% 1/2 current solution + 1/2 next solution
					%..........
					pop(n) = 0.5*(pop(n) + prod/(meta_loss + Kloss(n)));

				end % n

				% renormalize solution to statisfy condition sum(pop) == xn0.
				pop = (xn0/sum(pop))*pop;
				    
				% new timestep: calculate deviation
				%-----------------------------------
				xmaxdev=0.0;
				    
				for n = 1:Nuse
					if (pop(n) ~= 0.0) || (popold(n) ~= 0.0)
						dev=abs(pop(n) - popold(n))*2 / (pop(n) + popold(n));
					else
						dev = 0;
			 		end
	
					if dev > xmaxdev
						xmaxdev=dev;
					end
				end % n

				if xmaxdev < 1.e-6 %then we can exit the loop, we've reached an equilibrium
					disp(['Iterations: ' num2str(it)]);
					break
				end
			
			end % it

			if it == maxIterations
				disp ('max iterations reached. no solution found.');
				break
			end



			%keep desired levels from this solution
			if allCombos
				level_table(:, Te_i, ne_i, n0_i) = pop(start_level:end_level);
			else
				level_table(:, Te_i) = pop(start_level:end_level);
			end % allCombos

		end %n0_i
        
	end %ne_i

end %Te_i

% save solutions together with all solution parameters
solutions = struct('levels', level_table, 'Te', Te_list, 'ne', ne_list, 'n0', n0_list, 'ion_fraction', ion_fraction, 'plasma_radius', plasma_radius, 'combos', allCombos, 'start', start_level, 'end', end_level, 'Nuse', Nuse, 'useIonRecomb', useIonRecomb,'useRadDecay', useRadDecay, 'useThermal', useThermal, 'useDiff', useDiff, 'useMeta', useMeta);
