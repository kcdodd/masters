function [ integral ] = calc_integral(y_i, dx)
%CALC_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

    N = size(y_i,1);
    
    integral = zeros(N,1);
    integral(1) = 0;
    
    for i = 2:2:N
        
        %trapazoid_rule for middle element
        integral(i) = integral(i-1) + dx*0.5*(y_i(i-1) + y_i(i));
        
        if (i+1 <= N)
        
            %simpsons rule for end element
            integral(i+1) = integral(i-1) + 2*dx*(y_i(i-1) + 4*y_i(i) + y_i(i+1))/6;
            
        end
        
    end
end

