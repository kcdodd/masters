function [solution] = arI_solution(Te_in, ne_in, xn0, ArII_fraction, plasma_radius, rates, coronal)
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
useThermal = 0; % Use collisions with thermal atoms. 
useMeta = 0; % use losses between meta-stable levels and ground.
maxIterations = 1e6; % maximum number of iterations before giving up on solution

% neutral argon mass
m0 = 6.63e-26; 
%parameter(pi=3.1415926536)
pi=3.1415926536;

% INITIALIZATION - PARAMETERS
%============================
% Total number of level parents considered
Ntot=65;

%precomputed gain/loss rates
Kloss(Ntot)=zeros;
Kgain(Ntot, Ntot)=zeros;
Rgain(Ntot)=zeros;

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

%disp('Begining Precalculations.');

    
    %fprintf(['\tCurrent division: ' num2str(r_i) ' of ' num2str(numdivs+1) '\n']);

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
		Te = Te_in; % get electron temperature
        
        dTe = rates.Te(2) - rates.Te(1);
        Te_min = rates.Te(1);
        
        kioneIntegral(n) = get_interp(rates.ion(:,n), Te_min, dTe, Te);
        rrec1eIntegral(n) = get_interp(rates.recomb_2b(:,n), Te_min, dTe, Te);
        rrec2eIntegral(n) = get_interp(rates.recomb_3b(:,n), Te_min, dTe, Te);
		

		for m = 1:Ntot
            
			% only transition to other levels
			if m ~= n
                    		
				% m > n does not mean Eexc(m) > Eexc(n), so determine if this is excitation or 
				% de-excitation
				if (Eexc(m) > Eexc(n))
                    
                    kexceIntegral(n,m) = get_interp(rates.trans_e(:,n,m), Te_min, dTe, Te);


				else
					% this is a de-excitation transion
                    
                    kdeexeIntegral(n,m) = get_interp(rates.trans_e(:,n,m), Te_min, dTe, Te);

				end

			end
            
		end
        
	end% n
        

    % get current electron density
    nes = ne_in;

    % Define ion density. It can be different from ne if Nediv~=1
    ni=nes*ArII_fraction;

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
        Kloss(n) = useIonRecomb*(kione(n) + useThermal*kionth(n));

        % radiative, 3-body, and thermal recombination
        %- - - - - - - - - - - - - - - - - - - - - - - - - 
        Rgain(n) = useIonRecomb*(rrec1e(n) + rrec2e(n) + useThermal*rrecth(n));

        for m = 1:Nuse


            if m ~= n
                if (Eexc(m) > Eexc(n))

                    % excitation

                    % electron, ion, fast atom, therm.atom excitation to higher levels
                    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    % (escape factors to level 1 incorporated)

                    if n == 1
                        %Arad=A(n,m) * Esc(m);
                        Arad=A(n,m);
                    else
                        Arad=A(n,m);
                    end

                    Kloss(n) = Kloss(n) + kexce(n,m) + useThermal*kexcth(n,m);

                    Kgain(n, m) = kdeexe(m, n) + useThermal*kdeexth(m,n) + useRadDecay*Arad;
                else
                    %de-excitation

                    % elec,ion,fast atom,therm.atom deexcit, radiat.decay to lower levels
                    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    % (escape factors to level 1 incorporated)

                    if m == 1
                        %Arad=A(m,n) * Esc(n);
                        Arad=A(m,n);
                    else
                        Arad=A(m,n);
                    end

                    Kloss(n)= Kloss(n) + kdeexe(n,m) + useThermal*kdeexth(n,m) + useRadDecay*Arad;

                    Kgain(n, m) = kexce(m, n) + useThermal*kexcth(m,n);
                end
            end

        end % m

    end % n
    

%disp('Precalculations done. Begining iterations.');

% population:
npop = zeros(Nuse,1);


    % initial guess is all in ground state
    npop(1) = xn0;

if (coronal == 1)
        for n = 2:Nuse

            G = Kgain(n, 1)*npop(1);

            L = Kloss(n);

            npop(n) = G/L;
        end
else
    % solve next
    for it = 1:maxIterations

        npop_old = npop;

        for n = 2:Nuse

            G = Rgain(n);

            for m = 1:Nuse
                G = G + Kgain(n, m)*npop_old(m);
            end

            L = Kloss(n);

            npop(n) = G/L;
        end

        % new timestep: calculate deviation to see if solution is found
        %-----------------------------------
        xmaxdev=0.0;

        for n = 1:Nuse
            if (npop(n) ~= 0.0) || (npop_old(n) ~= 0.0)
                dev=abs(npop(n) - npop_old(n))*2 / (npop(n) + npop_old(n));
            else
                dev = 0;
            end

            if dev > xmaxdev
                xmaxdev=dev;
            end
        end % n


        if xmaxdev < 1.e-7 %then we can exit the loop, we've reached a solution

            break;
        end
    end

    if (it == maxIterations)
        fprintf(['\t\tNo Soltution Found (static): ' num2str(it) '\n']);

    end
    
end



% save solutions together with all solution parameters
solution = struct('npop', npop, 'Te', Te_in, 'ne', ne_in, 'ArII_fraction', ArII_fraction, 'Nuse', Nuse, 'useIonRecomb', useIonRecomb,'useRadDecay', useRadDecay, 'useThermal', useThermal, 'useMeta', useMeta);
