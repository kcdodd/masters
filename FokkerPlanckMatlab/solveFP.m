function [ evdf, Te_eff, velocities, K_exc, K_ion] = solveFP( Te, ne, N, rates, rates_energies, rates_delta_energy, varargin )
%SOLVEFP Summary of this function goes here
% Te, energies, and rates_energy in eV
% rates in 1/s
% ne in m^-3
% eedf in 1/J
% Te_eff in eV
% varargin = optionally set initial evdf to work from

    constants;
    
    max_velocity = 4*sqrt(2*const_e*Te/const_me);
    dv = max_velocity/N;
    
    velocities = (0:dv:max_velocity)';

    N = size(velocities, 1);
    num_rates = size(rates, 2);
    
    dv = (velocities(2) - velocities(1));
    dE = rates_energies(2) - rates_energies(1);

    if (length(varargin) == 1)
        
        evdf = varargin{1};
        
    else
        % start with maxwellian @ Te
        evdf = exp(-0.5*const_me*velocities.^2/(const_e*Te));


        A = calc_integral_total(evdf.*(4*pi*velocities.^2), dv);

        evdf = evdf/A;
    
    end
    
    rates_power = zeros(N,1);
    for i=1:N
        energy = 0.5*const_me*velocities(i)^2/const_e;
        
        for r = 1:num_rates
        
            rates_power(i) = rates_power(i) + get_interp(rates(:,r), rates_energies(1), dE, energy)*rates_delta_energy(r);
        end
    end
    
    max_iterations = 1;
    tolerance = 1E-4;
    
    for it = 1:max_iterations
        
        [evdf, Te_eff, K_exc, K_ion] = solveFP_0(evdf, Te, ne, velocities, rates, rates_energies, rates_delta_energy);
        
        if (2*abs(Te_eff - Te)/(Te_eff + Te) < tolerance)
            break;
        end
        
        % adjust R0 to try to match target temperature
        %R0 = R0*(1 + 2*(Te-Te_eff)/(Te_eff + Te));

    end
end

