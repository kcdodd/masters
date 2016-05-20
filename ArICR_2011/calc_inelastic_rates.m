function [ rates, rates_energy ] = calc_inelastic_rates(n0, energies, inelastic_corrections, varargin )
%CALC_INELASTIC_RATES Summary of this function goes here
% rates in 1/s
% energies in eV
% varargin = filename to export rates
    constants;
    
    %From Bogaerts model, Sciamma's code.
    
    Ntot=65;

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
    
    %use only rates from ground state
    
    %number of rates = 66 = 65 + 1 including ionization
    num_rates = 65;
    
    N = size(energies, 1);
    
    rates = zeros(N, num_rates);
    rates_energy = zeros(num_rates, 1);
    
	for m = 2:65
        rates_energy(m-1) = Eexc(m)-Eexc(1);
	end
    
    rates_energy(65) = Eion(1);
    
    for i = 1:N
    
        for m = 2:65
            
            if (energies(i) > rates_energy(m-1))

                % 1E-4 for cm^2 to m^2
                rates(i, m-1) = sqrt(2*const_e*energies(i)/const_me)*(1E-4)*Sexce(energies(i) ,1,m,Eexc(1),Eexc(m),aAf(1, m),aP(1, m),g(1),g(m));
            end
        end
        
        if (energies(i) > rates_energy(65))

            rates(i, 65) = sqrt(2*const_e*energies(i)/const_me)*(1E-4)*Sione(energies(i) ,1,Eion(1));
            
        end

    end
    
    rates = n0*rates;
    

    
    if (length(varargin) == 1)
        
        filename = varargin{1};
        
        fileID = fopen(filename, 'w');
        
        % output number of electron energy values
        fprintf(fileID, '%f\r\n', N);
        
        % output number of rates
        fprintf(fileID, '%f\r\n', num_rates);
        
        % output the electron energies used
        for i = 1:N

            fprintf(fileID, '%f\t', energies(i, 1));
            
        end
        
        fprintf(fileID, '%s\r\n', '');
        
        % output the delta energy for each rate
        for r = 1:num_rates

            fprintf(fileID, '%f\t', rates_energy(r, 1));
            
        end
        
        fprintf(fileID, '%s\r\n', '');
        
        % output all rates
        
        for r = 1:num_rates

            for i = 1:N

                fprintf(fileID, '%f\t', rates(i, r));
            end
            
            fprintf(fileID, '%s\r\n', '');
            
        end


        fclose(fileID);
    end
end








