function [sum] = sum_rates(velocities, evdf, rates, rates_energies, rates_delta_energy)
%SUM_RATES Summary of this function goes here
%   Detailed explanation goes here

    vmin = 9E4;

    constants;

    num_rates = size(rates, 2);
    
    N = size(velocities,1);
    
    loss_r = zeros(N,num_rates);
    
    dE = rates_energies(2)-rates_energies(1);
    dv = velocities(2)-velocities(1);
    
    i_min = floor(vmin/dv)+1;
    vmin = (i_min-1)*dv;
    
    for i = 1:N
        energy = 0.5*const_me*(velocities(i)^2)/const_e;
        
        for r = 1:num_rates
            
            loss_r(i,r) = (velocities(i)^2)*evdf(i)*get_interp(rates(:,r), rates_energies(1), dE, energy);
        end
    end
    
    loss = zeros(N,1);
    gain = zeros(N,1);
    
    
    for i = i_min:N
        energy = 0.5*const_me*(velocities(i)^2)/const_e;
        
        for r = 1:num_rates
            
            source_energy = energy + rates_delta_energy(r);
            source_velocity = sqrt(2*const_e*source_energy/const_me);
            
            loss(i) = loss(i) + loss_r(i,r);
            
            if (source_velocity <= velocities(N))

                gain(i) = gain(i) + get_interp(loss_r(:,r), velocities(1), dv, source_velocity);

            end

        end
    end
    
    for i = 1:(i_min-1)
        gain(i) = gain(i_min)*(velocities(i)/vmin)^2;
    end
    
    gain_t = calc_integral_total(gain, velocities(2)-velocities(1));
    loss_t = calc_integral_total(loss, velocities(2)-velocities(1));
    
    gain = gain*loss_t/gain_t;
    
    sum = (gain - loss)./(velocities.^2);
    
    for i = 1:(i_min-1)
        sum(i) = sum(i_min);
    end
    
end

