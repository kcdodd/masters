function [ out ] = Qdi(E, I, A, B, C, D)
%QDI Summary of this function goes here
%   Detailed explanation goes here
% E in eV
    n = size(E, 2);
    
    out = zeros(1, n);

    for i = 1:n
   
        
        if (E(i) < I)
            out(i) = 0;
        else

            u = E(i)/I;

            %[m^2]
            out(i) = (1E-4)*(1E-14)*(u*I^2)^(-1)*(A*(1-1/u) + B*(1-1/u)^2 + C*log(u) + D*log(u)/u);
        end

    end

end

