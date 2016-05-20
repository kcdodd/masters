function [ factor ] = calc_correction_factor(min_energy, f, f0, velocities)
%UNTITLED5 Summary of this function goes here
%   calculates rate correction factor

    constants;

    N = size(f, 1);

    int_f = calc_integral(4*pi*(velocities.^2).*f, velocities(2));
    int_f0 = calc_integral(4*pi*(velocities.^2).*f0, velocities(2));
    
    i_min = 0;
    
    %find velocity nearest min_energy
    for i = 1:N
        energy = 0.5*const_me*velocities(i)^2/const_e;
        
        if (energy > min_energy)
            i_min = i-1;
            break;
        end
    end
    
    if (i_min == 0)
       disp('did not find energy'); 
       
       return;
    end
    
    if (int_f(N) > 1)
        %disp(['not normal: ', num2str(int_f(N), '%f') ]); 
        %int_f
    end

    
    factor = (int_f(N) - int_f(i_min))/(1 - int_f0(i_min));
end

