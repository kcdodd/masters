function [exc, dexc, ion, rec, dec] = calc_rates(ne, Te) % ne in m-3, Te in Joules
n_max = 5;

osc_strengths = [
    0   2E-1    5E-2    2E-2  9E-3;
    0	0       6E-1    1E-1    1E-2;
    0   0       0       1       1E-1;
    0   0       0       0       1.3];

exc = zeros(n_max, n_max);
dexc = zeros(n_max, n_max);
ion = zeros(n_max, 1);
rec = zeros(n_max, 1);
dec = zeros(n_max, n_max);

me = 9.11E-31;
e = 1.602E-19;
a0 = 5.29E-11;
hbar = 1.05E-34;
h= 6.626E-34;
ttivity = 8.85E-12;
r0 = 2.82E-15;

vth = sqrt(Te/me);

mults = 1:1000;
dv = vth/10;

for i = 1:n_max
	E_i = -13.6*4*e/i^2;
    
	for j = (i+1):n_max
		E_j = -13.6*4*e/j^2;
		dE = E_j - E_i;
		v_0 = sqrt(2*dE/me);
		v = mults.*dv + v_0;
		G = ((2*pi*hbar^4*osc_strengths(i, j))/(me^3*a0^2*dE))*sqrt((2/pi)*(me/(Te))^3)*exp(-(me/(2*Te))*(v.^2));
		exc(i,j) = ne*trapz(v, G);
		dexc(j,i) = exc(i,j)*(i^2/j^2)*exp(dE/Te);

		dec(j, i) = osc_strengths(i, j)*2*r0*(dE/hbar)^2/(3E8);
	end

    
	v_0 = sqrt(2*(-E_i)/me);
	v = mults.*dv + v_0;
    
	G =	(e^4/(16*pi^2*ttivity^2*0.5*me))*(1/(-E_i)-1./(0.5*me*v.^2))*sqrt((2/pi)*(me/(Te))^3).*exp(-(me/(2*Te))*(v.^2));
	ion(i) = ne*trapz(v, G);

	
	rec(i) = ion(i)*ne*(i^2)*(h^2/(me*Te))^(3/2)*exp(E_i/Te);
end
