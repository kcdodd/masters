function [Te_best, n0_best, emission_best, emission_actual, lines, levels_best, levels_actual, R2] = computeRadialProfile_2(wavelengths, spectra, ne_m, Te_l, n0_l, level_table)

Nspecs = size(spectra, 2);

Te_best = zeros(Nspecs,1);
n0_best = zeros(Nspecs,1);


for i = 1:Nspecs
	[Te_best(i), n0_best(i), emission_best(:,i), emission_actual(:,i), lines, levels_best(:,i), levels_actual(:,i), R2(i)] = computeBestFitCR_4(wavelengths, spectra(:,i),  ne_m(i), Te_l, n0_l, level_table);
end