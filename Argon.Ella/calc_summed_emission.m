function [sums_exp, sums_calc, lines] = calc_summed_emission(wavelengths, exp_emission, calced_lines)
lineData = [
    696.5431 3;
    706.7218 3;
    738.3980 3;
    751 3;
    763.5106 3;
    772.4 2;
    794.8176 3;
    801 3;
    811 4;
    826.4522 3;
    841.5 3;
    852.1442 3;
    912.2967 3;
    922.4499 3;
    965.7786 3];

sums_exp = zeros(size(lineData, 1), 1);
lines = zeros(size(lineData, 1), 1);

sums_calc = [
    calced_lines(1);
    calced_lines(2);
    calced_lines(3);
    calced_lines(4) + calced_lines(5);
    calced_lines(6);
    calced_lines(7) + calced_lines(8);
    calced_lines(9);
    calced_lines(10) + calced_lines(11);
    calced_lines(12) + calced_lines(13);
    calced_lines(14);
    calced_lines(15) + calced_lines(16);
    calced_lines(17);
    calced_lines(18);
    calced_lines(19);
    calced_lines(20)];
    

for i = 1:size(lineData, 1)
    lines(i) = lineData(i, 1);
    
    for w = 1:size(wavelengths, 1)
        if (wavelengths(w) > lineData(i, 1) - lineData(i, 2) && wavelengths(w) < lineData(i, 1) + lineData(i, 2))
            sums_exp(i) = sums_exp(i) + exp_emission(w);
        end
    end
    
end