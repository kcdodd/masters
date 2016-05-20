function [ evdf , Te_eff, K_exc, K_ion] = solveFP_0(evdf_0, Te, ne, velocities, rates, rates_energies, rates_delta_energy)
%SOLVEFP_0 Summary of this function goes here
% Te, energies, and rates_energy in eV
% eedf_0 in 1/J
% R0 in eV/s
% rates in 1/s
% ne in m^-3
% eedf in 1/J
% Te_eff in eV

    constants;

    N = size(velocities, 1);
    
    num_rates = size(rates, 2);
    
    dv = (velocities(2) - velocities(1));
    
    %dE = rates_energies(2) - rates_energies(1);
    
    evdf = evdf_0;
    
    max_iterations = 1;
    tolerance = 1E-8;
    dt = 1e-10;
    
    collision_factor = ne*(gamma_ee(Te, ne)/(const_me^2));

    heating_factor = 1;
    
    for it = 1:max_iterations
        D = collision_factor*calc_D(velocities, evdf);
        R = collision_factor*calc_R(velocities, evdf);
        D_derivative = calc_derivative(D, dv);
        R_derivative = calc_derivative(R, dv);
        
        %K_exc = sum_rates(velocities, evdf, rates(:,1:(num_rates-1)), rates_energies, rates_delta_energy(1:(num_rates-1)));
        K_exc = calc_excitations_1(velocities, evdf, rates(:,1:(num_rates-1)), rates_energies, rates_delta_energy(1:(num_rates-1)));
        
        K_ion = calc_ionization_1(velocities, evdf, rates(:,num_rates), rates_energies, rates_delta_energy(num_rates));
        
        K = K_exc + K_ion;

        
        evdf_derivative = calc_derivative(evdf, dv);
        evdf_derivative_2 = calc_derivative_2(evdf, dv);
        
        M = evdf_derivative*2./velocities;
        M(1)=M(2);
        
        calc_integral_total(K.*(4*pi*velocities.^4), dv)
        
        D0 = -calc_integral_total(K.*(4*pi*velocities.^4), dv)/calc_integral_total((M + evdf_derivative_2).*(4*pi*velocities.^4), dv);
        
        M = (evdf_derivative.*(2*(D + heating_factor*D0)))./velocities;
        M(1) = M(2);

        M2 = (2*R)./velocities;
        M2(1) = M2(2);

        k = (M + (D_derivative - R).*evdf_derivative + evdf_derivative_2.*(D + heating_factor*D0) - (M2 + R_derivative).*evdf + K);

        evdf_next = evdf + dt*k;

        evdf_next = max(evdf_next, 0);


       
        A = calc_integral_total(evdf_next.*(4*pi*velocities.^2), dv);
        
        evdf_next = evdf_next/A;
        
        avg_dev = calc_integral_total((4*pi*velocities.^2).*abs(evdf_next-evdf), dv);
        
        evdf = evdf_next;
        
        Te_eff = (1/3)*const_me*calc_integral_total(evdf.*(4*pi*velocities.^4), dv)/const_e;
        
        heating_factor = Te/Te_eff;
        
        disp(['End iteration ', num2str(it), ', Te = ', num2str(Te_eff), ', avg_dev = ', num2str(avg_dev)]);
        
        if (avg_dev < tolerance)
            %converged
            disp('converged');
            break;
        end
        
    end
    
    %calculate new effective temperature
    
end

