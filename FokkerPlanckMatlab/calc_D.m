function [ D ] = calc_D(velocities, evdf)
%CALC_D Summary of this function goes here
%   Detailed explanation goes here

    N = size(velocities, 1);
    
    dv = (velocities(2)-velocities(1));
    
    D = zeros(N,1);
    
    a = evdf.*velocities;
    
    for i = 1:N
        velocity = velocities(i);
        
        f = a.*G(velocities/velocity);
        
        D(i) = 2*pi*calc_integral_total(f, dv);
        
    end
    
end

