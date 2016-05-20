function [ out ] = G(beta)
%G Summary of this function goes here
%   Detailed explanation goes here

    N=size(beta,1);
    
    out = zeros(N,1);
    
    for i = 1:N
        
        if (beta(i) < 1)
            out(i) = beta(i)^3*4/3;
        else
            out(i) = 4/3;
        end
    
    end
end

