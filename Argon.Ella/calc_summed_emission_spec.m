function [sums_spec, lines] = calc_summed_emission_spec(wavelengths, emission)
lineData = [
    696.5431 3; %1
    706.7218 3; %2
    714.7042 3; %3
    727.2936 3; %4
    738.3980 3; %5
    751 3; %6
    763.5106 3; %7
    772.4 2; %8
    794.8176 3; %9
    801 3; %10
    811 4; %11
    826.4522 3; %12
    841.5 4; %13
    852.1442 3; %14
    866.7943 3; %15
    912.2967 3; %16
    922.4499 3; %17
    935.4220 3; %18
    965.7786 3; %19
    978.4503 3]; %20

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
