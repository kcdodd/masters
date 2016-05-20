function [ power, profile ] = RadiativeCooling( Te, ne, n0, r , V)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    coeffs = [-3.7573E-20, 2.8398E-17, -7.4088E-15, 2.2873E-12, -1.3614E-10, 2.4647E-9];
    
    power = 0;
    weights = 0;
    
    n_r = max(size(r,1), size(r,2));
    profile = zeros (1, n_r);
    
    for r_i = 1:n_r
        
        if (Te(r_i) > 2)
            
        
            p_i = 0;

            for i=0:5
                p_i = p_i + coeffs(i+1)*(Te(r_i)/1000)^i;
            end

            profile(r_i) = p_i*(1E-7)*(1E-6);
            power = power + p_i*(1E-7)*(1E-6)*ne(r_i)*n0(r_i)*r(r_i);
            weights = weights + r(r_i);
        end
    end
    
    
    power = V*power/weights;
end

