function int = calcLampIntegration(wavelengths, temp, power, time, ratio)

wavelengths_nm = wavelengths*(1E-9);
% dlambda is 1nm
int = ratio.*(1E-9)*time*power*15*(6.626E-34*3E8)^3./(pi^4*(1.38E-23*temp)^4*(wavelengths_nm.^4).*(exp(6.626E-34*3E8./(wavelengths_nm*1.38E-23*temp))-1));