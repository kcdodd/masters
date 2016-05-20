function [rates, r_divs] = ionrates(Te_r_in, ne_r_in, r_in, rmin, rmax, numdivs, element, outfile)

r_divs = zeros(1, numdivs+1);

rates = zeros(1, numdivs+1);

dr = (rmax-rmin)/numdivs;

% compute radii
for i = 1:(numdivs+1)
    r_divs(i) = dr*(i-1)+rmin;
end

% interp ne_r
ne_r = interp1(r_in, ne_r_in, r_divs, 'linear');

%interp Te_r
Te_r = interp1(r_in, Te_r_in, r_divs, 'linear');

found = 0;


if (strcmp(element, 'H+0'))
    C = Cdi_Hydrogen_I(Te_r);
    found = 1;
end

if (strcmp(element, 'He+0'))
    C = Cdi_Helium_I(Te_r);
    found = 1;
end

if (strcmp(element, 'Ne+0'))
    C = Cdi_Neon_I(Te_r);
    found = 1;
end

if (strcmp(element, 'Ar+0'))
    C = Cdi_Argon_I(Te_r);
    found = 1;
end

if (found == 0)
    disp(['Element ' element ' not found.']);
    return;
end

%size(C)
%size(ne_r)
rates = ne_r.*C;

fileID = fopen(outfile, 'w');

fprintf(fileID, '%u\r\n', (numdivs+1));

for nr = 1:(numdivs+1)
    fprintf(fileID, '%f\t', r_divs(nr));
end

fprintf(fileID, '%s\r\n', '');

for nr = 1:(numdivs+1)
    fprintf(fileID, '%f\t', Te_r(nr));
end

fprintf(fileID, '%s\r\n', '');

for nr = 1:(numdivs+1)
    fprintf(fileID, '%E\t', ne_r(nr));
end

fprintf(fileID, '%s\r\n', '');

for nr = 1:(numdivs+1)
    
    fprintf(fileID, '%E\r\n', rates(nr));
end


fclose(fileID);

