function Srec2e_A = Srec2e(E_A,n,Eion,gn,gin, xkT)
% electron 3b recombination (from ionization; inverse process)

% KCD vectorized E and return.

% !!! Srec2ne = sigma/n_e !!!
% ! 3/2 kT=mean energy of all electrons, assume: Ee=10 eV
% ! fac32=(h**2/(2pi*M_e))**3/2 : in eV**3/2 (cfr kT) cm**3 (cfr n_e)

N_E = size(E_A, 2);

Srec2e_A = zeros(1, N_E);

for E_i = 1:N_E
    
    E = E_A(E_i);

    %parameter(Ee=10.0,fac32=3.313e-22)
    %Ee=10;
    fac32=3.313e-22;
    Srec2e=0.0;

    %if(E.gt.0.0)then
    if E > 0.0
        %xkT=0.667 * Ee;
        fac32b=fac32 / xkT^1.5;
        E2=E + Eion;

        Srec2e=(gn/(2*gin)) * fac32b * (E2/E) * Sione(E2,n,Eion);
    end
    
    Srec2e_A(E_i) = Srec2e;
end % E_i

%return
%end
