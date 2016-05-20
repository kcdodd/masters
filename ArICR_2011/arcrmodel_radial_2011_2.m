function [solution] = arcrmodel_radial_2011_2(Te_r_in, ne_r_in, r_in, xn0, T0, ion_fraction, rmin, rmax, numdivs)
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
Nuse=35; % how many levels to include: 1 up to max is 65.
useIonRecomb = 1; % use ionization and recombination. =1 to use process, 0 to not use it.
useRadDecay = 1; % use spontaneous radiative decay
useThermal = 0; % Use collisions with thermal atoms. 
useMeta = 0; % use losses between meta-stable levels and ground.
useQuadGK = true; % set to false if quadgk is not supported
maxIterations = 1e6; % maximum number of iterations before giving up on solution

% neutral argon mass
m0 = 6.63e-26; 
%parameter(pi=3.1415926536)
pi=3.1415926536;

% INITIALIZATION - PARAMETERS
%============================
% Total number of level parents considered
Ntot=65;

% Population variables
dr = (rmax - rmin)/numdivs;
v_th = 100*sqrt(2*T0*1.6e-19/m0);
dt = dr/v_th;
r_divs(numdivs+1, 1) = zeros;
ne_r(numdivs+1, 1) = zeros;
Te_r(numdivs+1, 1) = zeros;
plasma_radius = (rmax-rmin)/2;

% compute radii
for i = 1:(numdivs+1)
    r_divs(i) = dr*(i-1)+rmin;
end

% interp ne_r
ne_r = interp1(r_in, ne_r_in, r_divs, 'linear');

%interp Te_r
Te_r = interp1(r_in, Te_r_in, r_divs, 'linear');

%precomputed gain/loss rates
Kloss(numdivs+1, Ntot)=zeros;
Kloss_nexc(numdivs+1, Ntot)=zeros;
Kloss_ion(numdivs+1, Ntot)=zeros;
Kgain(numdivs+1, Ntot, Ntot)=zeros;
Rgain_ion(numdivs+1, Ntot)=zeros;

% Metastable diffusion, 2- and 3-body recombination coefficient rates
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

% Escape factor
Esc(Ntot,1)=zeros;

% metastable variables
bbi(4,1)=zeros;
beta=zeros;
gamma=zeros;

%=================================================================
% DATA INPUT
%============

% xtgas is the ambient gas temperature in K 
xtgas=348;
pres = 1;

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

disp('Begining Precalculations.');

for r_i = 1:(numdivs+1)
    
    fprintf(['\tCurrent division: ' num2str(r_i) ' of ' num2str(numdivs+1) '\n']);

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
		Te = Te_r(r_i); % get electron temperature

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
        

    % get current electron density
    nes = ne_r(r_i);

    % Define ion density. It can be different from ne if Nediv~=1
    ni=nes*ion_fraction;

    % multiply all integrals by factors of ni and ne. now k**** has units of cm^3/s, and r**** units of 1/s.
    kione = nes*kioneIntegral;
    rrec1e = ni*nes*rrec1eIntegral;
    rrec2e = ni*nes*nes*rrec2eIntegral;
    kexce = nes*kexceIntegral;
    kdeexe = nes*kdeexeIntegral;

    kionth=zeros(Ntot,1);
    rrecth=zeros(Ntot,1);
    kexcth=zeros(Ntot,Ntot);
    kdeexth=zeros(Ntot,Ntot);

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
    for n = 1:Nuse

        % electron, ion, fast atom, therm.atom ionization
        %- - - - - - - - - - - - - - - - - - - - - - - - - 
        Kloss_ion(r_i, n) = useIonRecomb*(kione(n) + useThermal*kionth(n));

        % radiative, 3-body, and thermal recombination
        %- - - - - - - - - - - - - - - - - - - - - - - - - 
        Rgain_ion(r_i, n) = useIonRecomb*(rrec1e(n) + rrec2e(n) + useThermal*rrecth(n));

        for m = 1:Nuse


            if m ~= n
                if (Eexc(m) > Eexc(n))

                    % excitation

                    % electron, ion, fast atom, therm.atom excitation to higher levels
                    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    % (escape factors to level 1 incorporated)

                    if n == 1
                        Arad=A(n,m) * Esc(m);
                    else
                        Arad=A(n,m);
                    end
                    
                    if (m == 1 || m == 2 || m == 4)
                        Kloss_nexc(r_i, n) = Kloss_nexc(r_i, n) + kexce(n,m) + useThermal*kexcth(n,m);
                    end

                    Kloss(r_i, n) = Kloss(n) + kexce(n,m) + useThermal*kexcth(n,m);

                    Kgain(r_i, n, m) = kdeexe(m, n) + useThermal*kdeexth(m,n) + useRadDecay*Arad;
                else
                    %de-excitation

                    % elec,ion,fast atom,therm.atom deexcit, radiat.decay to lower levels
                    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    % (escape factors to level 1 incorporated)

                    if m == 1
                        Arad=A(m,n) * Esc(n);
                    else
                        Arad=A(m,n);
                    end
                    
                    if (m == 1 || m == 2 || m == 4)
                        Kloss_nexc(r_i, n) = Kloss_nexc(r_i, n) + kdeexe(n,m) + useThermal*kdeexth(n,m) + useRadDecay*Arad;
                    end

                    Kloss(r_i, n)= Kloss(r_i, n) + kdeexe(n,m) + useThermal*kdeexth(n,m) + useRadDecay*Arad;

                    Kgain(r_i, n, m) = kexce(m, n) + useThermal*kexcth(m,n);
                end
            end

        end % m

    end % n
    
end %r_i

disp('Precalculations done. Begining iterations.');


% n > 4
nexc = zeros(2, numdivs+1);

% 1=ground, 2=(n=2), 3=(n=4), 4=(n>4)
dndt = zeros(2, numdivs+1, 4);

% population: 1=particles moving +R, 2=-R
npop = zeros(2, numdivs+1, Nuse);

% boundary conditions
npop(1, 1, 1) = 0.5*xn0;
npop(2, numdivs+1, 1) = 0.5*xn0;
dndt(1, 1, 1) = -npop(1, 1, 1)*Kloss_ion(1, 1);
dndt(2, numdivs+1, 1) =  -npop(2, numdivs+1, 1)*Kloss_ion(numdivs+1, 1);

for r_i = 2:(numdivs+1)
    % find solution for current division
    
    fprintf(['\tCurrent division: ' num2str(r_i) ' of ' num2str(numdivs+1) '\n']);
    
    rcur = [r_i, numdivs+2-r_i];
    rprev = [r_i-1, numdivs+3-r_i];
    
    for p = 1:2
        
        % initial guess is from last r
        npop(p, rcur(p), :) = npop(p, rprev(p), :);

        % solve positive npop
        for it = 1:maxIterations

            npop_old = npop(p, rcur(p), :);

            % compute ground state
            G=0;

            for m = 2:Nuse
                G = G + Kgain(rcur(p), 1, m)*npop(p, rcur(p), m);
            end

            L = Kloss_ion(rcur(p), 1) + Kloss(rcur(p), 1);

            npop(p, rcur(p), 1) = (0.5*dt*(dndt(p, rprev(p), 1) + G) + npop(p, rprev(p), 1))/(1 + 0.5*dt*L);
            dndt(p, rcur(p), 1) = G - npop(p, rcur(p), 1)*L;

            % compute meta-stable levels
            %+++++++++++++
            for n = [2 4]

                % additional losses for meta-stable levels (meta-meta collisions + meta-ground collisions)
                %...............................................
                if n == 2
                    pop2=npop_old(4);
                end

                if n == 4
                    pop2=npop_old(2);
                end

                meta_loss = 2*xkmet*npop_old(n) + xkmet*pop2 + npop_old(1)*k2b(n) + npop_old(1)^2*k3b(n);


                G = 0;

                for m = 1:Nuse
                    G = G + Kgain(rcur(p), n, m)*npop_old(m);
                end

                L = Kloss(rcur(p), n) + meta_loss + Kloss_ion(rcur(p), n);

                if n == 2
                    npop(p, rcur(p), 2) = (0.5*dt*(dndt(p, rprev(p), 2) + G) + npop(p, rprev(p), 2))/(1 + 0.5*dt*L);

                    if (npop(p, rcur(p), 2) < 0)
                        npop(p, rcur(p), 2) = 0;
                    end

                    dndt(p, rcur(p), 2) = G - L*npop(p, rcur(p), 2);
                end

                if n == 4
                    npop(p, rcur(p), 4) = (0.5*dt*(dndt(p, rprev(p), 3) + G) + npop(p, rprev(p), 4))/(1 + 0.5*dt*L);

                    if (npop(p, rcur(p), 4) < 0)
                        npop(p, rcur(p), 4) = 0;
                    end

                    dndt(p, rcur(p), 3) = G - L*npop(p, rcur(p), 4);
                end

            end

            
            G_exc = 0;
            L_exc = 0;

            % compute excited levels
            %+++++++++++++
            for n = [3 5:Nuse]

                G_exc = G_exc + Kgain(rcur(p), n, 1)*npop_old(1) + Kgain(rcur(p), n, 2)*npop_old(2) + Kgain(rcur(p), n, 4)*npop_old(4);
                L_exc = L_exc + npop(p, rcur(p), n)*(Kloss_nexc(rcur(p), n) + Kloss_ion(rcur(p), n));

            end
            
            nexc_tot = sum(npop(p, rcur(p), [3 5:Nuse]));

            if (nexc_tot ~= 0)
                L_exc = L_exc/nexc_tot;
            else
                L_exc = 0;
            end

            nexc(p, rcur(p)) = (0.5*dt*(dndt(p, rprev(p), 4) + G_exc) + nexc(p, rprev(p)))/(1 + 0.5*dt*L_exc);

            if (nexc(p, rcur(p)) < 0)
                nexc(p, rcur(p)) = 0;
            end
            
            for n = [3 5:Nuse]

                G = 0;

                for m = 1:Nuse
                    G = G + Kgain(rcur(p),n,m)*npop_old(m);
                end

                L = Kloss(rcur(p), n) + Kloss_ion(rcur(p), n);

                if (L ~= 0)
                    npop(p, rcur(p), n) = G/L;
                else
                    fprintf(['\t\tL=0 (+): ' num2str(n) '\n']);
                end

            end
            
            nexc_tot = sum(npop(p, rcur(p), [3 5:Nuse]));

            if (nexc_tot ~= 0)
                npop(p, rcur(p), [3 5:Nuse]) = npop(p, rcur(p), [3 5:Nuse])*nexc(p, rcur(p))/nexc_tot;
            else
                fprintf('\t\tnexc_tot = 0\n');
            end

            dndt(p, rcur(p), 4) = G_exc - L_exc*nexc(p, rcur(p));

            % new timestep: calculate deviation to see if solution is found
            %-----------------------------------
            xmaxdev=0.0;

            for n = 1:Nuse
                if (npop(p, rcur(p), n) ~= 0.0) || (npop_old(n) ~= 0.0)
                    dev=abs(npop(p, rcur(p), n) - npop_old(n))*2 / (npop(p, rcur(p), n) + npop_old(n));
                else
                    dev = 0;
                end

                if dev > xmaxdev
                    xmaxdev=dev;
                end
            end % n


            if xmaxdev < 1.e-7 %then we can exit the loop, we've reached an equilibrium
                fprintf(['\t\tIterations (+): ' num2str(it) '\n']);
                break
            end
        end
    end
    
end

vr = zeros(numdivs+1);
npop_tot = zeros(numdivs+1, Nuse);

for r_i = 1:(numdivs+1)
    npop_tot(r_i, :) = npop(1, r_i, :) + npop(2, r_i, :);
    
    sum_pos = sum(npop(1, r_i, :));
    sum_neg = sum(npop(2, r_i, :));
    
    vr(r_i) = (sum_pos - sum_neg)/(sum_pos + sum_neg);
end

vr = v_th*vr;

% save solutions together with all solution parameters
solution = struct('population', npop_tot, 'velocity', vr, 'npop', npop, 'nexc', nexc, 'dndt', dndt, 'Te', Te_r, 'ne', ne_r, 'r', r_divs, 'ion_fraction', ion_fraction, 'rmin', rmin, 'rmax', rmax, 'numdivs', numdivs, 'Nuse', Nuse, 'useIonRecomb', useIonRecomb,'useRadDecay', useRadDecay, 'useThermal', useThermal, 'useMeta', useMeta);
