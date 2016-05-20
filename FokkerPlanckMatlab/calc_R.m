function [ R ] = calc_R(velocities, evdf)
%CALC_R Summary of this function goes here
%   Detailed explanation goes here

    N = size(velocities, 1);
    
    dv = (velocities(2)-velocities(1));
    
    R = zeros(N,1);
    
    evdf_derivative = calc_derivative(evdf, dv);
    
    for i=1:N
        velocity = velocities(i);
        
        f = evdf_derivative.*G(velocities/velocity);

        R(i) = 2*pi*velocity*calc_integral_total(f, dv);
    end
    


end

