function [ K ] = calc_excitations_1(  velocities, evdf, excitation_rates, rate_energies, excitation_energies  )
%CALC_EXCITATIONS_1 Summary of this function goes here
%   Detailed explanation goes here


    constants;

    N = size(velocities, 1);
    
    num_energies = size(rate_energies,1);
    num_rates = size(excitation_rates,2);
    
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
    
    L = zeros(num_energies,1);
    
    for e = 1:num_energies
           
        L(e) = eedf(e)*sum(excitation_rates(e, :));
        
    end
    
    %calc gain rate
    G = zeros(num_energies,1);
    
    
    for e = 1:num_energies
        for r = 1:num_rates
            source_energy = rate_energies(e) + excitation_energies(r);
            
            if (source_energy < rate_energies(num_energies))
           
                G(e) = G(e) + get_interp(eedf, rate_energies(1), dE, source_energy)*get_interp(excitation_rates(:, r), rate_energies(1), dE, source_energy);
            end
        end
    end
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

