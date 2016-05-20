function [ out ] = gamma_ee(Te, ne)
%GAMMA_EE Summary of this function goes here
%   Detailed explanation goes here

    constants;

    coulomb_logarithm = 23.5 - 0.5*log(ne/1E6) + (5/4)*log(Te) - sqrt(1E-5 + (log(Te)-2)^2/16);
    

    
    out = const_e^4*coulomb_logarithm/(8*pi*const_epsilon_0^2);
end

