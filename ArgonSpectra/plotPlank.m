function [plankCurve, ratio] = plotPlank(wavelengths, curve, zero, temp)

wavelengths_nm = wavelengths*(1E-9);

plankCurve = curve(zero)*(wavelengths_nm(zero)^4)*(exp(6.626E-34*3E8/(wavelengths_nm(zero)*1.38E-23*temp))-1)./((wavelengths_nm.^4).*(exp(6.626E-34*3E8./(wavelengths_nm*1.38E-23*temp))-1));

plot(wavelengths, curve, wavelengths, plankCurve);

ratio = curve./plankCurve;