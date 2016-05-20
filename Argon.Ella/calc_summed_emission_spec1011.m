function [sums_spec, lines] = calc_summed_emission_spec1011(wavelengths, emission)
lineData = [
    751 3];

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
