function [ h ] = Fit_1(R, Ra, Rb, A, Rp, gamma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xp = (Rp - Ra)/(Rb - Ra);

x = (R - Ra)/(Rb - Ra);

y = x.^2*(0.5-xp)/(xp^2 - xp) + x*(1-(0.5-xp)/(xp^2 - xp));

h = A*((gamma/2)*(1-cos(pi*y)) + (1-gamma/2)*sin(pi*y).^3);
end

