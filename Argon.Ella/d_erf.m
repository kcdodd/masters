function d_erf = d_erf(E)
% error function in code did not compile, added 9/22/03 by AMK
% 07-04-18 by EMS replaced "derf" by "d_erf" since DERF is a built
% fortran function
	
d_erf=0.0;

if E < 0
	d_erf=-gammp(.5,E^2);
else
	d_erf=gammp(.5,E^2);
end
        
%	return
%	end
