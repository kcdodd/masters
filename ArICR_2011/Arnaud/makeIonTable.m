function [ ion_table ] = makeIonTable( input_args )
%MAKEIONTABLE Summary of this function goes here
%   Detailed explanation goes here
Te_l = [0.5 1 2 5 10 20];%eV

ion_table = struct('Te', Te_l, 'C1', Cdi_2(Te_l, 1), 'C2', Cdi_2(Te_l, 2));

end

