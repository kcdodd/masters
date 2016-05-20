function [xn0,ArI_pop,ArI_Gain, ArI_Loss] = arcrmodel_new3(nes, Te, xn0, Nediv, xrad)
%==========================================================
% Ar I Collisional radiative model (1D)
% received by Ella Sciamma(EMS)from Amy Keesee
% modified October 2007 by EMS 
% --> convertion from fortran to matlab
% --> took off radial dependency
% 
% get ne and Te from "ne_Te_filename"
% "factor" is used to calculate the neutral density from ne
% --> if factor=0.25, n0=ne*0.25 => degree of ionization
%       ne/(ne+n0)=80%
% "xrad" is the radius of the plasma column
% --> used in the calculation of the metastable diffusion
% "press" is the fill pressure
% --> 0.05 mTorr for VX-100
% --> 1 mTorr for Helicon old and new and Helimak
%==========================================================

% INITIALIZATION - PARAMETERS
%============================
% Total number of level parents considered
Ntot=65;

% Which things to include in model
Nuse=65; %how many levels to include: 1 up to max is 65.
useIonRecomb = 1; % use ionization and recombination
useRadDecay = 1; % use spontaneous radiative decay
useThermal = 1; % Use collisions with thermal atoms. 1 to use process, 0 to not use it.
useDiff = 1; % use diffusion of meta-stable levels
useMeta = 1; % use losses between meta-stable levels.

% Number of energy point considered for EEDF
NE=739;

% Population variable
pop(Nuse,1)=zeros;

% Metastable diffusion, 2- and 3-body recombination coefficient rates
D(4)=zeros;
k2b(4)=zeros;
k3b(4)=zeros;

% electron density "nes" and ion density "ni"
%nes=zeros;
ni=zeros;

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
% EEDF
fe(NE+1,1)=zeros;

% rate coefficients
kione(Ntot,1)=zeros;
rrec1e(Ntot,1)=zeros;
rrec2e(Ntot,1)=zeros;
kexce(Ntot,Ntot,1)=zeros;
kdeexe(Ntot,Ntot,1)=zeros;
kioni(Ntot,1)=zeros;
rreci(Ntot,1)=zeros;
kexci(Ntot,Ntot,1)=zeros;
kdeexi(Ntot,Ntot,1)=zeros;
kiona(Ntot,1)=zeros;
rreca(Ntot,1)=zeros;
kexca(Ntot,Ntot,1)=zeros;
kdeexa(Ntot,Ntot,1)=zeros;
kionth(Ntot,1)=zeros;
rrecth(Ntot,1)=zeros;
kexcth(Ntot,Ntot)=zeros;
kdeexth(Ntot,Ntot)=zeros;
Kloss(Ntot,1)=zeros;

% metastable variables
bbi(4,1)=zeros;
beta=zeros;
gamma=zeros;

% neutral density
%xn0=zeros;

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
        min=n;
        max=Ntot-n;
    else
        min=min + Ntot - (n-1);
        max=max + Ntot - n;
    end
    aAf(n,n+1:Ntot)=levels2m(min:max,3)';
    aP(n,n+1:Ntot)=levels2m(min:max,4)';
    A(n,n+1:Ntot)=levels3m(min:max,3)';
end


    % EEDF.DAT
    % - - - - -
%Electron energy
Ee(1:11)=0:0.001:0.01;
Ee(12:20)=0.02:0.01:0.1;
Ee(21:29)=0.2:0.1:1;
Ee(30:43)=1.2:0.2:3.8;
Ee(44:NE+1)=4:1:700; 

% Get electron density "nes" and EEDF "fe" "calcul_EEDF.m"
%[fe]= calcul_EEDF_2(nes,Te,Ee);
[fe]= calcul_EEDF_3(Te,Ee);
fe = nes*fe;


% Define ion density. It can be different from ne if Nediv~=1
ni=nes/Nediv;


%4. calculation of ground state from ne
%---------------------------------------
%xn0= nes*factor;

%====================================================================
% CALCULATION OF COLLISION 'RATES' (r: rate (in cm-3 s-1), k: per pop=1 (in s-1)
%=================================
% ELECTRONS
%----------
for n = 1:Ntot
    %disp(['n= ' num2str(n)]) % to follow the evolution of the calculation on matlab
    
    for iEe = 1:NE+1
        kione(n)=kione(n) + fe(iEe) * Sione(Ee(iEe),n,Eion(n));
        
        rrec1e(n)=rrec1e(n) + fe(iEe) * Srec1e(Ee(iEe),n,Eion(n),gam(n),g(n),gi(n)) * fac(n) * ni;

        rrec2e(n)=rrec2e(n) + fe(iEe) * Srec2e(Ee(iEe),n,Eion(n),g(n),gi(n), Te) * nes * fac(n) * ni;
        
        for m=n+1:Ntot % for n=Ntot, the code will just not go through this loop
            kexce(n,m)=kexce(n,m) + fe(iEe) * Sexce(Ee(iEe),n,m,Eexc(n),Eexc(m),aAf(n,m),aP(n,m),g(n),g(m));
        end
        
        for l = 1:n-1
            kdeexe(n,l)=kdeexe(n,l) + fe(iEe) * Sdeexe(Ee(iEe),n,l,Eexc(n),Eexc(l),aAf(l,n),aP(l,n),g(n),g(l),g0(n));
        end       
    end
end

% THERMALIZED ATOMS (E=0.03 eV -> v=3.81e4 cm/s) (k: per pop=1, in s-1)
%-------------------
Eth=0.03; % thermal energy
vth=3.81e4; % thermal velocity

for n = 1:Ntot
    kionth(n)=Sionth(Eion(n)) * vth * xn0;
    
    rrecth(n)=Sreci(Eth,n,Eion(n),g(n),gi(n)) * vth * xn0 * nes * fac(n) * ni; 
    
    for m = n+1:Ntot
        kexcth(n,m)=Sexcth(n,m,Eexc(n),Eexc(m),gi(n),gi(m)) * vth * xn0;
    end
    
    for l = 1:n-1
        kdeexth(n,l)=Sdeexth(n,l,Eexc(n),Eexc(l),g(n),g(l),gi(n),gi(l),g0(n)) * vth * xn0;
    end
end

%==========================================================================
% CALCULATION OF THE ESCAPE FACTORS
%==================================
for n = 2:Ntot
    if A(1,n) ~= 0.0
        %xkR=2.1e-17/(Eexc(n)^3) * g(n)/(xtgas^0.5) * A(1,n) * xn0 * xrad;
        %xa=A(1,n) * (1 + 3.225e-14/(Eexc(n)^3) * g(n) * xn0) * 4.839e-9/Eexc(n)/(xtgas^0.5);
        %Td=1/(xkR * (pi*log(xkR))^0.5);
        %Tc=(xa / (pi^0.5 * xkR))^0.5;
        %Tcd=2 * xa / (pi * (log(xkR))^0.5);
        
        %Esc(n)=1.9 * Td * exp(-pi * Tcd^2 / (4.0*Tc^2)) + 1.3 * Tc * d_erf(pi^0.5 * Tcd/(2.0*Tc));
        %Esc(n)=real(Esc(n));

	xkR=(2.1e-17)* g(n) * A(1,n) * xn0 * xrad/((Eexc(n)^3)*sqrt(xtgas));
	xa=A(1,n) * (1 + (3.225e-14)* g(n) * xn0/(Eexc(n)^3))*(4.839e-9)/(Eexc(n)*sqrt(xtgas));
	Tc=sqrt(xa/(sqrt(pi)*xkR));

	if (xkR > 1.0)
		
		Td=1/(xkR*sqrt(pi*log(xkR)));
		Tcd=2*xa / (pi*sqrt(log(xkR)));

		Esc(n) = 1.9*Td*exp(-pi*(Tcd^2)/(4*Tc^2))+1.3*Tc*erf(sqrt(pi)*Tcd/(2*Tc));
	else
		Esc(n) = 1.3*Tc;
	end

	if Esc(n) > 1.0
		Esc(n) = 1.0;
	end
    else
        Esc(n)=1.0;
    end
end

%==================================================================
% CALCULATION OF THE LEVEL POPULATIONS
%======================================
% Initial values: all populations=0, pop(ground n=1)=xn0
%----------------------------------------------------------
pop(1)=xn0;
pop(2:Nuse)=zeros;

losses = zeros(Nuse, Nuse+1);

% ct terms (indep.of time) for all levels: loss by :
%----------------------------------------------------
for n = 2:Nuse
    Kloss(n)=0.0;
    
    % electron, ion, fast atom, therm.atom excitation to higher levels
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for m = n+1:Nuse

            
        Kloss(n)=Kloss(n) + kexce(n,m) + useThermal*(kexci(n,m) + kexca(n,m) + kexcth(n,m));
        losses(n, m) = kexce(n,m) + useThermal*(kexci(n,m) + kexca(n,m) + kexcth(n,m));
        
    end
    
    % elec,ion,fast atom,therm.atom deexcit, radiat.decay to lower levels
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % (escape factors to level 1 incorporated)
    for l = 1:n-1
       
        if l == 1
            Arad=A(l,n) * Esc(n);
        end
        
        if l ~= 1
            Arad=A(l,n);
        end
        
        losses(n, l) = useRadDecay*Arad + kdeexe(n,l) + useThermal*(kdeexi(n,l) + kdeexa(n,l) + kdeexth(n,l));
        
        Kloss(n)=Kloss(n) + kdeexe(n,l) + useThermal*(kdeexi(n,l) + kdeexa(n,l) + kdeexth(n,l)) + useRadDecay*Arad;
        
    end
    
    % electron, ion, fast atom, therm.atom ionization
    %- - - - - - - - - - - - - - - - - - - - - - - - - 
    Kloss(n)=Kloss(n) + useIonRecomb*(kione(n) + kioni(n) + kiona(n) + useThermal*kionth(n));
    
    losses(n, Nuse+1) = useIonRecomb*(kione(n) + kioni(n) + kiona(n) + useThermal*kionth(n));
end


% ct.terms (indep.of time) for metast: 2,4
%-----------------------------------------

for n = 2:2:4
    Kloss(n)=Kloss(n) + useMeta*(xn0*k2b(n) + xn0^2*k3b(n));
    %losses(n, Nuse+1) = losses(n, Nuse+1) + xn0*k2b(n) + xn0^2*k3b(n);
    
    bbi(n)=1/dt + useDiff*(D(n)/xrad * (2/xrad)) + Kloss(n);
    %bbi(n)=1/dt + Kloss(n);
end


%Build Rate Matrix
%------------------

%AMat = zeros(Ntot,Ntot);
%BVec = zeros(Ntot, 1);

%for n = 1:Ntot
%    BVec(n) =  -(rrec1e(n) + rrec2e(n) + rreci(n) + rreca(n) + rrecth(n));
    
%    for i = 1:Ntot
        
%        if i < n
%            AMat(n,i) = kexce(i,n) + kexci(i,n) + kexca(i,n) + kexcth(i,n);
%        elseif (i == n && i > 1)

            % electron, ion, fast atom, therm.atom ionization
            %- - - - - - - - - - - - - - - - - - - - - - - - - 
%             AMat(n,n) = - (sum(losses(n,:) + kione(n) + kioni(n) + kiona(n) + kionth(n)));
%        elseif i > n
%            AMat(n,i) = kdeexe(i,n) + kdeexi(i,n) + kdeexa(i,n) + kdeexth(i,n) + A(n,i);
%        end
            
%    end
%end




% Each timestep (pop(1->n-1) known at t, pop(n+1->Ntot) known at t-1)
%--------------------------------------------------------------------
mult=10;
t=0.0;

ArI_Gain = zeros(Nuse, Nuse+1);
ArI_Loss = zeros(Nuse, Nuse+1);

for it = 1:10000000
    t=t+dt;
    
    % old values
    %- - - - - - -
    popold(2:Nuse)=pop(2:Nuse); %vectorized loop
    
    % each level
    %- - - - - - -
    
        % metastable levels: 2,4
        %++++++++++++++++++++++++

    for n = 2:2:4
        
        % additional loss (met-met coll; incorpor.in bi)
        %...............................................
        if n == 2
            pop2=pop(n+2);
        end
        
        if n == 4
            pop2=pop(n-2);
        end
        
        bi=bbi(n) + useMeta*(2*xkmet*pop(n) + xkmet*pop2);
        %bi=bbi(n);
        
        prod=0.0;
        
        % prod: excit.from lower levels
        %..............................
        for l = 1:n-1
            prod=prod + (kexce(l,n) + useThermal*(kexci(l,n) + kexca(l,n) + kexcth(l,n))) * pop(l);
            ArI_Gain(n, l) = (kexce(l,n) + useThermal*(kexci(l,n) + kexca(l,n) + kexcth(l,n))) * pop(l);
            ArI_Loss(n,l) =  losses(n, l)*pop(n);
            
        end
        
        % prod: deexcit.+rad.decay from higher levels
        %............................................
        for m=n+1:Nuse
            prod=prod + (kdeexe(m,n) + useThermal*(kdeexi(m,n) + kdeexa(m,n) + kdeexth(m,n)) + useRadDecay*A(n,m))*pop(m);
            ArI_Gain(n, m) = (kdeexe(m,n) + useThermal*(kdeexi(m,n) + kdeexa(m,n) + kdeexth(m,n)) + useRadDecay*A(n,m))*pop(m);
            ArI_Loss(n,m) = losses(n, m)*pop(n);

        end
        
        % prod: recomb  -> total prod
        %.............................
        prod=prod + useIonRecomb*(rrec1e(n) + rrec2e(n) + rreci(n) + rreca(n) + rrecth(n));
        
        ArI_Gain(n, Nuse+1) = useIonRecomb*(rrec1e(n) + rrec2e(n) + rreci(n) + rreca(n) + rrecth(n));
        ArI_Loss(n, Nuse+1) = bi*pop(n);
        
        di=popold(n)/dt + prod;
        
        beta=bi;
        gamma=di;
        pop(n)=gamma/beta;
        
    end

        % other levels
        %+++++++++++++
    for n = 3:Nuse
        if n ~= 4 %we took care of n=2 and n=4 already
            
            % production: sum over all levels
            %..................................
            prod=0.0;
            
            % prod: excit.from lower levels
            %...............................
            for l=1:n-1
                %if (l ~= 2 && l ~= 4)
                    
                prod=prod + (kexce(l,n) + useThermal*(kexci(l,n) + kexca(l,n) + kexcth(l,n))) * pop(l);
                ArI_Gain(n, l) = (kexce(l,n) + useThermal*(kexci(l,n) + kexca(l,n) + kexcth(l,n))) * pop(l);
                ArI_Loss(n,l) = losses(n, l)*pop(n);
                
                %end
            end
            
            % prod: deexcit.+rad.decay from higher levels
            %.............................................
            for m = n+1:Nuse
                
                %if (m ~= 4)
                    
                prod=prod + (kdeexe(m,n) + useThermal*(kdeexi(m,n) + kdeexa(m,n) + kdeexth(m,n)) + useRadDecay*A(n,m))*pop(m);
                ArI_Gain(n,m) = (kdeexe(m,n) + useThermal*(kdeexi(m,n) + kdeexa(m,n) + kdeexth(m,n)) + useRadDecay*A(n,m))*pop(m);
                ArI_Loss(n,m) = losses(n, m)*pop(n);
                
                
                %end
            end
            
            % prod: recomb  -> total prod
            %.............................
            prod=prod + useIonRecomb*(rrec1e(n) + rrec2e(n) + rreci(n) + rreca(n) + rrecth(n));
            
            ArI_Gain(n, Nuse+1) = useIonRecomb*(rrec1e(n) + rrec2e(n) + rreci(n) + rreca(n) + rrecth(n));
            ArI_Loss(n, Nuse+1) = losses(n, Nuse+1)*pop(n);
            
            % equation
            %..........
            pop(n)=(pop(n) + dt*prod) / (1 + Kloss(n)*dt);
        end
    end
    
    % new timestep: calculate deviation
    %-----------------------------------
    xmaxdev=0.0;
    
    for n = 2:Nuse
        if (pop(n) ~= 0.0) || (popold(n) ~= 0.0)
            dev=abs(pop(n) - popold(n))*2 / (pop(n) + popold(n));
        else
            dev = 0;
        end
        
        if xmaxdev >= dev
            xmaxdev=xmaxdev;
        else
            xmaxdev=dev;
        end
    end
    
    if it == mult
        mult=mult+10;
    end
    
    if xmaxdev < 1.e-6 %then we can exit the loop, we've reached an equilibrium
        break
    end
end





if xmaxdev < 1.e-6
    %disp('t is')
    %disp(t)
    %disp('xmaxdev is less than 1e-6 => equilibrium')
    %disp('xmaxdev is ')
    %disp(xmaxdev);
end

n2=1:Nuse;

ArI_pop=pop;




%ArI_pop = linsolve(AMat, BVec);

%for n = 1:Ntot
%    for i = 1:Ntot
%        if (n ~= i)
%            ArI_Contributions(n,i) = AMat(n,i)*ArI_pop(i);
%        end
%    end
%end

    


%normalize each level's contributions
for n=1:Nuse
    totalRate = sum(ArI_Gain(n,:));
    
    if (totalRate ~= 0.0)
        ArI_Gain(n,:) = ArI_Gain(n,:)/totalRate;
    end
    
    totalRate = sum(ArI_Loss(n,:));
    
    if (totalRate ~= 0.0)
        ArI_Loss(n,:) = ArI_Loss(n,:)/totalRate;
    end
end
