function [ out ] = Cdi_Hydrogen_I(kT)
%CDI Summary of this function goes here
%   Detailed explanation goes here
% 1 = Ar+0
% 2 = Ar+1

    n = size(kT, 2);
    
    out = zeros(1,n);

    params = [
        13.6, 22.8, -12.0, 1.9, -22.6];
    
    for i = 1:n
        if (kT(i) == 0)
            out(i) = 0;
        else

            x_1 = params(1, 1)/kT(i);


            F_1 = F(x_1, params(1, 2), params(1, 3), params(1, 4), params(1, 5));

            %[m^3 s^-1]
            out(i) = (1E-6)*(6.69E-7)*(kT(i)^(-3/2))*(exp(-x_1)*F_1/x_1);
        end
    end
end

