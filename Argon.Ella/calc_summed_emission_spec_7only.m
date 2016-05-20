function [sums_spec, lines] = calc_summed_emission_spec_7only(wavelengths, emission)
lineData = [
    763.5106 3;
    801 3;
    811 4;
    866.7943 3;
    922.4499 3;
    935.4220 3;
    978.4503 3];

sums_spec = zeros(size(lineData, 1), 1);
lines = zeros(size(lineData, 1), 1);

    

for i = 1:size(lineData, 1)
    lines(i) = lineData(i, 1);
    
    for w = 1:size(wavelengths, 1)
        if (wavelengths(w) > lineData(i, 1) - lineData(i, 2) && wavelengths(w) < lineData(i, 1) + lineData(i, 2))
            sums_spec(i) = sums_spec(i) + emission(w);
        end
    end
    
end
