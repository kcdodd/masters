function [ flux ] = calc_flux( R0, Te, ne, energies, eedf)
%CALC_FLUX Summary of this function goes here
%   Detailed explanation goes here

    constants;

    N = size(energies, 1);
    flux = zeros(N,1);
    
    eedf_derivative = calc_derivative(eedf, energies);
    
    for i=1:N
        energy = energies(i);
        flux(i) = calc_D(i, energies, eedf, Te, ne)*(2*const_me*energy*eedf_derivative(i) + const_me*eedf(i)) - (calc_R(i, energies, eedf, Te, ne)*sqrt(2*energy/const_me) + R0)*eedf(i);
    end
end

