function [ out ] = Cdi_2(Te, Ar_state )
%CDI_2 Summary of this function goes here
%   Detailed explanation goes here
% kT in eV
%Ar_state: 1=Ar+0, 2=Ar+1

    me = 9.11E-31; %electron mass is kg

        params = [
            15.8, 171.1, -78.0, 3.8, -169.0;
            27.6, 147.0, -97.4, 3.2, -107.7;
            40.9, 122.8, -81.4, 2.6, -90.0];
    
    E_0 = 0.001;
    E_1 = 0.01;
    E_2 = 0.1;
    E_3 = 1;
    E_4 = 4;
    E_5 = 20;
    E_6 = 100;
    E_7 = 1000;
    
    n = size(Te, 2);
    
    out = zeros(n,1);

    for i = 1:n
        
        % define the v*f part of sigma*v*f to compute
        % <sigma*v>: f=(3D isotropic maxwell-boltzmann). with v in m/s, Te and Ee in eV.
        vf = @(Ee) 1.6E-19*sqrt((2/pi)*(me/(1.6E-19*Te(i)))^3).*(2*1.6E-19*Ee/(me*me)).*exp(-Ee/Te(i));

        % define electron impact ionization integrand
        kioneU = @(Ee) vf(Ee).*Qdi(Ee, params(Ar_state, 1), params(Ar_state, 2), params(Ar_state, 3), params(Ar_state, 4), params(Ar_state, 5)); % define sgima*v*f product

        out(i) = quadgk(kioneU, params(Ar_state, 1), params(Ar_state, 1)+E_1)+quadgk(kioneU, params(Ar_state, 1)+E_1, params(Ar_state, 1)+E_2)+quadgk(kioneU, params(Ar_state, 1)+E_2, params(Ar_state, 1)+E_3)+quadgk(kioneU, params(Ar_state, 1)+E_3, params(Ar_state, 1)+E_4)+quadgk(kioneU, params(Ar_state, 1)+E_4, params(Ar_state, 1)+E_5)+quadgk(kioneU, params(Ar_state, 1)+E_5, params(Ar_state, 1)+E_6)+quadgk(kioneU, params(Ar_state, 1)+E_6, params(Ar_state, 1)+E_7);
    end
end

