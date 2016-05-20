function [ K ] = calc_ionization_1( velocities, evdf, ionization_rate, rate_energies, ionization_energy )
%CALC_IONIZATION Summary of this function goes here
%   Detailed explanation goes here

    constants;

    N = size(velocities, 1);
    
    num_energies = size(rate_energies,1);
    
    dv = (velocities(2) - velocities(1));
    dE = rate_energies(2) - rate_energies(1);
    
    eedf = zeros(num_energies,1);

    %convert to eedf from evdf
    for e = 1:num_energies
        v = sqrt(2*rate_energies(e)*const_e/const_me);
        
        if (v <= velocities(N))
        
            eedf(e) = (4*pi*const_e*v/const_me)*get_interp(evdf, velocities(1), dv, v);
        else
            eedf(e) = 0;
        end
        
    end
    
    % make sure it's still normalized in terms of numerical integral
    A = calc_integral_total(eedf, dE);
    eedf = eedf/A;
    
    
    %calc gain rate
    G = zeros(num_energies,1);
    
    
    for e = 1:num_energies
        
        gain_density = eedf.*ionization_rate./(rate_energies - ionization_energy);
        
        for e2 = 1:num_energies
            if (rate_energies(e2) < ionization_energy + rate_energies(e))
                
                gain_density(e2) = 0;
            else
                break;
            end
        end
        
        % factor of 2 is because of contribution equally from freed electron and incident electron
        G(e) = 2*calc_integral_total(gain_density, dE);
    end
    
    % loss rate is only the total rate at each energy, notice gain is 2x
    % higher than loss, which gives electron production.
    
    L = eedf.*ionization_rate;
    
    % add artificial "diffusive" loss L_0 to compensate for Sum(Gain) =
    % 2xSum(L) such that Sum(Gain) = Sum(L + L_0)
    
    
    L_0 = calc_integral_total(L, dE);
    
    L = eedf.*(ionization_rate + L_0);
    
    % now convert back to evdf version and ensure balance once again
    
    G_v = zeros(N,1);
    L_v = zeros(N,1);
    
    for i = 2:N
        energy = 0.5*const_me*velocities(i)^2/const_e;
        
        G_v(i) = (const_me/(4*pi*const_e*velocities(i)))*get_interp(G, rate_energies(1), dE, energy);
        L_v(i) = (const_me/(4*pi*const_e*velocities(i)))*get_interp(L, rate_energies(1), dE, energy);
    end
    
    G_v(1) = G_v(2);
    L_v(1) = L_v(2);
    
    G_v = smooth(G_v);
    
    L_0_v = calc_integral_total(4*pi*(velocities.^2).*L_v, dv);
    G_0_v = calc_integral_total(4*pi*(velocities.^2).*G_v, dv);
    
    G_v = G_v*L_0_v/G_0_v;
    
    K = G_v - L_v;
    
end

