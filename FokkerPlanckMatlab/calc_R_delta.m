function [ R ] = calc_R_delta( velocities, v0 )
%CALC_R_DELTA Summary of this function goes here
%   Detailed explanation goes here

    R = -2*pi*G_prime(v0./velocities);
end

