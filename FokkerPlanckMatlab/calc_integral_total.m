function [ integral ] = calc_integral_total(y_i, dx)
%CALC_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

    N = size(y_i,1);
    
    integral = 0;
    
    for i = 2:2:N
        
        if (i+1 <= N)
        
            %simpsons rule
            integral = integral + 2*dx*(y_i(i-1) + 4*y_i(i) + y_i(i+1))/6;
            
        end
        
    end
    
    if (mod(N,2) == 0)
        %have to do the last point by trapazoid if even number of points
        integral = integral + 0.5*dx*(y_i(N-1) + y_i(N));
    end
end

