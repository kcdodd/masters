function [n0_best, emission_best, emission_actual, lines, levels_best, levels_actual, R2] = computeRadialProfile(wavelengths, spectra, n0_l, level_tables)

Nspecs = size(spectra, 2);

Te_best = zeros(Nspecs,1);
n0_best = zeros(Nspecs,1);


for i = 1:Nspecs
	[n0_best(i), emission_best(:,i), emission_actual(:,i), lines, levels_best(:,i), levels_actual(:,i), R2(i)] = computeBestFitCR_3(wavelengths, spectra(:,i), n0_l, level_tables(:,:,i));
end
