function [ out ] = G_prime( beta )
%G_PRIME Summary of this function goes here
%   Detailed explanation goes here

    N=size(beta,1);
    
    out = zeros(N,1);
    
    for i = 1:N
        
        if (beta(i) < 1)
            out(i) = (beta(i)^2)*4;
        else
            out(i) = 0;
        end
    
    end
end

