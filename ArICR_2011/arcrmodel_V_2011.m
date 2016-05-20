function [solution] = arcrmodel_V_2011(Te_r_in, ne_r_in, r_in, xn0, T0, ion_fraction, rmin, rmax, z_max, numdivs, rates)
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

% Population variables
dr = (rmax - rmin)/numdivs;
dz = zmax/numdivs;
v_th = 100*sqrt(2*T0*1.6e-19/m0);
dt_min = 2e-7;
dt = dr/v_th;
plasma_radius = (rmax-rmin)/2;

numrin = max(size(r_in,1),size(r_in,2));

%precomputed gain/loss rates
Kloss(numrin, Ntot)=zeros;
Kgain(numrin, Ntot, Ntot)=zeros;
Rgain(numrin, Ntot)=zeros;

% solution grid
n_solution = zeros(numdivs+1, numdivs+1, Nuse);

r_grid = zeros(numdivs+1);
z_grid = zeros(numdivs+1);

for i = 1:numdivs+1
    r_grid(i) = (i-1)*dr + rmin;
    z_grid(i) = (i-1)*dz;
end

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

for r_i = 1:numrin
    
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
		Te = Te_r_in(r_i); % get electron temperature
        
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
    nes = ne_r_in(r_i);

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
        Kloss(r_i, n) = useIonRecomb*(kione(n) + useThermal*kionth(n));

        % radiative, 3-body, and thermal recombination
        %- - - - - - - - - - - - - - - - - - - - - - - - - 
        Rgain(r_i, n) = useIonRecomb*(rrec1e(n) + rrec2e(n) + useThermal*rrecth(n));

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

                    Kloss(r_i, n) = Kloss(r_i, n) + kexce(n,m) + useThermal*kexcth(n,m);

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

                    Kloss(r_i, n)= Kloss(r_i, n) + kdeexe(n,m) + useThermal*kdeexth(n,m) + useRadDecay*Arad;

                    Kgain(r_i, n, m) = kexce(m, n) + useThermal*kexcth(m,n);
                end
            end

        end % m

    end % n
    
end %r_i

disp('Precalculations done. Begining iterations.');

% define the directions from which to sample incoming particle paths
directions = [
	1, 0, 0;
    1, -1, 0;
    0, -1, 0;
    -1, -1, 0;
    -1, 0, 0;
	1, 0, 1;
    1, -1, 1;
    0, -1, 1;
    -1, -1, 1;
    -1, 0, 1;
	1, 0, -1;
    1, -1, -1;
    0, -1, -1;
    -1, -1, -1;
    -1, 0, -1;
    0, 0, 1;
    0, 0, -1];

num_d = size(directions,1);

velocities = [v_th/4, v_th/2, v_th*3/4, v_th, v_th*5/4, v_th*3/2, v_th*7/4, v_th*2];

num_v = size(velocities, 2);

weights = zeros(size(velocities));

% calulate velocity weights
for v_i = 1:num_v
    weights(v_i) = exp(-velocities(v_i)^2/v_th^2);
end

% normalize weights
weights(:) = weights(:)/sum(weights(:));

% normalize directions
for d_i = 1:num_d
    directions(d_i, :) = directions(d_i, :)/sum(directions(d_i, :));
end

for d_i = 1:num_d
    
    for v_i = 1:num_v
        % create array of radius values
        
        dl = velocities(v_i)*dt_min;
    
    end
end
    
dndt = zeros(2, numdivs+1, Nuse);
dndt_ion = zeros(2, numdivs+1, Nuse);

% population: 1=particles moving +R, 2=-R
npop = zeros(2, numdivs+1, Nuse);
static = zeros(numdivs+1, Nuse);

% boundary conditions
npop(1, 1, 1) = 0.5*xn0;
npop(2, numdivs+1, 1) = 0.5*xn0;
dndt(1, 1, 1) = Rgain(1, 1) - npop(1, 1, 1)*Kloss(1, 1);
dndt_ion(1, 1, 1) = Rgain(1, 1) - npop(1, 1, 1)*Kloss(1, 1);
dndt(2, numdivs+1, 1) =  Rgain(numdivs+1, 1) - npop(2, numdivs+1, 1)*Kloss(numdivs+1, 1);
dndt_ion(2, numdivs+1, 1) = Rgain(numdivs+1, 1) - npop(2, numdivs+1, 1)*Kloss(numdivs+1, 1);

for n = 2:Nuse

    dndt(1, 1, n) = Rgain(1, n) + Kgain(1, n, 1)*npop(1, 1, 1);
    dndt_ion(1, 1, n) = Rgain(1, n);
    
    dndt(2, numdivs+1, n) = Rgain(numdivs+1, n) + Kgain(numdivs+1, n, 1)*npop(2, numdivs+1, 1);
    dndt_ion(2, numdivs+1, n) = Rgain(numdivs+1, n);
end

failed = false;

for r_i = 2:(numdivs+1)
    
    % find solution for current division
    
    %fprintf(['\tCurrent division: ' num2str(r_i) ' of ' num2str(numdivs+1) '\n']);
    
    rcur = [r_i, numdivs+2-r_i];
    rprev = [r_i-1, numdivs+3-r_i];
    
    for p = 1:2
        
        % initial guess is from last r
        npop(p, rcur(p), :) = npop(p, rprev(p), :);

        % solve next npop
        for it = 1:maxIterations

            npop_old = npop(p, rcur(p), :);

            for n = 1:Nuse

                G = Rgain(rcur(p), n);

                for m = 1:Nuse
                    G = G + Kgain(rcur(p), n, m)*npop_old(m);
                end

                L = Kloss(rcur(p), n);
                
                npop(p, rcur(p), n) = (0.5*dt*(dndt(p, rprev(p), n) + G) + npop(p, rprev(p), n))/(1 + 0.5*dt*L);
                npop(p, rcur(p), n) = max(npop(p, rcur(p), n), 0);
                
                dndt(p, rcur(p), n) = G - npop(p, rcur(p), n)*L;
            end
            
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


            if xmaxdev < 1.e-7 %then we can exit the loop, we've reached a solution
                if (p == 1)
                    %fprintf(['\t\tIterations (+): ' num2str(it) '\n']);
                else
                    %fprintf(['\t\tIterations (-): ' num2str(it) '\n']);
                end
                
                break
            end
        end
        
        if (it == maxIterations)
            if (p == 1)
                fprintf(['\t\tNo Soltution Found (+): ' num2str(it) '\n']);
            else
                fprintf(['\t\tNo Soltution Found (-): ' num2str(it) '\n']);
            end
            
            failed = true;

            break;
        end
    end
    
    if (failed == true)
        break;
    end

    % initial guess is all in ground state
    static(r_i, 1) = xn0;

    % solve next
    for it = 1:maxIterations

        npop_old = static(r_i, :);

        for n = 2:Nuse

            G = Rgain(r_i, n);

            for m = 1:Nuse
                G = G + Kgain(r_i, n, m)*npop_old(m);
            end

            L = Kloss(r_i, n);

            static(r_i, n) = G/L;
        end

        % new timestep: calculate deviation to see if solution is found
        %-----------------------------------
        xmaxdev=0.0;

        for n = 1:Nuse
            if (static(r_i, n) ~= 0.0) || (npop_old(n) ~= 0.0)
                dev=abs(static(r_i, n) - npop_old(n))*2 / (static(r_i, n) + npop_old(n));
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

        break;
    end
    
end

vr = zeros(numdivs+1,1);
npop_tot = zeros(numdivs+1, Nuse);
ntot = zeros(numdivs+1, 1);

for r_i = 1:(numdivs+1)
    npop_tot(r_i, :) = npop(1, r_i, :) + npop(2, r_i, :);
    
    sum_pos = sum(npop(1, r_i, :));
    sum_neg = sum(npop(2, r_i, :));
    ntot(r_i) = sum_pos + sum_neg;
    
    vr(r_i) = (sum_pos - sum_neg)/(sum_pos + sum_neg);
end

vr = v_th*vr;

% save solutions together with all solution parameters
solution = struct('density', ntot, 'population', npop_tot, 'velocity', vr, 'npop', npop, 'dndt', dndt, 'static', static, 'Te', Te_r_in, 'ne', ne_r_in, 'r', r_in, 'ion_fraction', ion_fraction, 'rmin', rmin, 'rmax', rmax, 'numdivs', numdivs, 'dt', dt, 'Nuse', Nuse, 'useIonRecomb', useIonRecomb,'useRadDecay', useRadDecay, 'useThermal', useThermal, 'useMeta', useMeta);