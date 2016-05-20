function [times, pop] = sim_pop(ne_0, Te_0, max_time, num_steps)

pop = zeros(5, num_steps);
dt = max_time/(num_steps-1);

times = dt*(0:(num_steps -1));

HeIII = ne_0/2;
ne = ne_0;

for s = 2:num_steps
	[exc, dexc, ion, rec, dec] = calc_rates(ne, Te_0);

	for i = 1:5

		dndt = 0;

		for j = 1:(i-1)
			dndt = dndt + pop(j, s-1)*exc(j, i) - pop(i, s-1)*(dexc(i, j) + dec(i, j));
		end

		
		for j = (i+1):5
			dndt = dndt + pop(j, s-1)*dexc(j, i) - pop(i, s-1)*exc(i, j);
		end

		dnIIIdt = -(HeIII*rec(i) - pop(i, s-1)*ion(i));

		dndt = dndt - dnIIIdt;

		pop(i, s) = pop(i, s-1) + dt*dndt;

		

		HeIII = HeIII + dt*dnIIIdt;

		ne = ne + dt*dnIIIdt;

		

	end
	
end

