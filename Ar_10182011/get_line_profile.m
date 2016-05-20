function profile = get_line_profile(wl, spec, center, width)

N_spec = size(spec,1);
N_wl = size(spec,2);

profile = zeros(N_spec, 1);

for s = 1:N_spec

	for i = 1:N_wl

		if (wl(i) > center - width && wl(i) < center + width)
			profile(s) = profile(s) + spec(s, i);
		end
	end
end
