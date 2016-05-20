function [ D ] = calc_D_delta(velocities,  v0)
%CALC_D_DELTA Summary of this function goes here
%   Detailed explanation goes here

	D = 2*pi*v0*G(v0./velocities);

end

