function [Te_best, n0_best, sums_calc_best, actual_emission, lines, n_levels_best, n_levels_actual, R2] = computeBestFitCR_4(wavelengths, spectral_emission, ne_i, Te_l, n0_l, level_table)


% wavelength A(i,j) level ratio
lineData = [
    696.5431 6.39e+06 4 1; % 1
    706.7218 3.80e+06 3 5/8; % 2
    714.7042 6.25e+5 3 5/8; % 3
    727.2936 1.83e6 4 1; % 4
    738.3980 8.47e+06 3 5/8; % 5
    750.3869 4.45e+07 6 1; % 6
    751.4652 4.02e+07 5 1; % 7
    763.5106 2.45e+07 2 5/20; % 8
    772.3761 5.18e+06 2 3/20; % 9
    772.4207 1.17e+07 4 1; % 10
    794.8176 1.86e+07 3 3/8; % 11
    800.6157 4.90e+06 2 5/20; % 12
    801.4786 9.28e+06 2 5/20; % 13
    810.3693 2.5e+07 2 3/20; % 14
    811.5311 3.31e+07 2 7/20; % 15
    826.4522 1.53e+07 4 1; % 16
    840.8210 2.23e+07 3 5/8; % 17
    842.4648 2.15e+07 2 5/20; % 18
    852.1442 1.39e+07 3 3/8; % 19
    866.7943 2.43e6 2 3/20; % 20
    912.2967 1.89e+7 1 1; % 21
    922.4499 5.03e+6 2 5/20; % 22
    935.4220 1.06e6 2 3/20; % 23
    965.7786 5.43e6 1 1; % 24
    978.4503 1.47e6 2 5/20]; % 25

lineGroupings = [
    696.5431 3 3 1 1; %1
    706.7218 3 3 2 2; %2
    714.7042 3 3 3 3; %3
    727.2936 3 3 4 4; %4
    738.3980 3 3 5 5; %5
    751 3 3 6 7; %6
    763.5106 3 3 8 8; %7
    772.4 2 2 9 10; %8
    794.8176 3 3 11 11; %9
    801 3 3 12 13; %10
    811 4 4 14 15; %11
    826.4522 3 3 16 16; %12
    841.5 4 4 17 18; %13
    852.1442 3 3 19 19; %14
    866.7943 3 3 20 20; %15
    912.2967 3 3 21 21; %16
    922.4499 3 3 22 22; %17
    935.4220 3 3 23 23; %18
    965.7786 3 3 24 24; %19
    978.4503 3 3 25 25]; %20
    
	actual_emission = zeros(size(lineGroupings, 1), 1);
	lines = zeros(size(lineGroupings, 1), 1);

	for i = 1:size(lineGroupings, 1)
	    lines(i) = lineGroupings(i, 1);
	    
	    for w = 1:size(wavelengths, 1)
		if (wavelengths(w) > lineGroupings(i, 1) - lineGroupings(i, 2) && wavelengths(w) < lineGroupings(i, 1) + lineGroupings(i, 3))
		    actual_emission(i) = actual_emission(i) + spectral_emission(w);
		end
	    end
	    
	end

n_levels_actual = zeros(6,1);

for n = 1:6
	n_sum = 0;
	norm = 0;

	for i = 1:size(lineGroupings, 1)
		allSame = true;
		denominator = 0;

		for j = lineGroupings(i,4):lineGroupings(i,5)
			if(lineData(j, 3) == n)
				denominator = denominator + lineData(j,2)*lineData(j,4);
			else
				allSame = false;
			end
		end %j

		if allSame
			n_sum = n_sum + actual_emission(i)^2/denominator;
			norm = norm + actual_emission(i);
		end
	end %i

	if (norm > 0)
		n_levels_actual(n) = n_sum/norm;
	end
end %n
    

    n0_best_i = 0;
	Te_best_i = 0;
    SSerr_best = Inf;
    
%    n_levels_actual = zeros(4,1);
    
    %calculate population levels from spectral lines
    % l = 6-9
%    n_levels_actual(1) = (actual_emission(16)^2/(lineData(21,2)*lineData(21,4)) + actual_emission(19)^2/(lineData(24,2)*lineData(24,4)))/(actual_emission(16) + actual_emission(19));
%    n_levels_actual(2) = (actual_emission(7)^2/(lineData(8,2)*lineData(8,4)) + actual_emission(10)^2/(lineData(12,2)*lineData(12,4) + lineData(13,2)*lineData(13,4)) + actual_emission(11)^2/(lineData(14,2)*lineData(14,4) + lineData(15,2)*lineData(15,4)) + actual_emission(15)^2/(lineData(20,2)*lineData(20,4)) + actual_emission(17)^2/(lineData(22,2)*lineData(22,4)) + actual_emission(18)^2/(lineData(23,2)*lineData(23,4)) + actual_emission(20)^2/(lineData(25,2)*lineData(25,4)))/(actual_emission(7)+actual_emission(10)+actual_emission(11)+actual_emission(15)+actual_emission(17)+actual_emission(18)+actual_emission(20));
%    n_levels_actual(3) = (actual_emission(2)^2/(lineData(2,2)*lineData(2,4)) + actual_emission(3)^2/(lineData(3,2)*lineData(3,4)) + actual_emission(5)^2/(lineData(5,2)*lineData(5,4)) + actual_emission(9)^2/(lineData(11,2)*lineData(11,4)) + actual_emission(14)^2/(lineData(19,2)*lineData(19,4)))/(actual_emission(2)+actual_emission(3)+actual_emission(5)+actual_emission(9)+actual_emission(14));
%    n_levels_actual(4) = (actual_emission(1)^2/(lineData(1,2)*lineData(1,4)) + actual_emission(4)^2/(lineData(4,2)*lineData(4,4)) + actual_emission(12)^2/(lineData(16,2)*lineData(16,4)))/(actual_emission(1)+actual_emission(4)+actual_emission(12));
    


        for Te_i = 1:size(Te_l, 2)
            for n0_i = 1:size(n0_l, 2)
                

                SSerr = 0;
                
                for n=1:6
			if (n_levels_actual(n) > 0)
                    		SSerr = SSerr + (level_table(n, Te_i, ne_i, n0_i) - n_levels_actual(n))^2;
			end
                end

                    
               
                if (SSerr < SSerr_best)
                    SSerr_best = SSerr;
                    n0_best_i = n0_i;
			Te_best_i = Te_i;
                end
                   

            end %n0
	end
        


Te_best = Te_l(Te_best_i);
    n0_best = n0_l(n0_best_i);
    
    mean = sum(n_levels_actual)/4;
    
    SStot = 0;
    
    n_levels_best = zeros(4,1);
    
    for i=1:4
        SStot = SStot + (n_levels_actual(i) - mean)^2;
        n_levels_best(i) = level_table(i, n0_best_i);
    end
    
    R2 = 1 - SSerr_best/SStot;
    
                    %Calculate lines from solution
                calced_lines = zeros(25, 1);
                for i = 1:25
                    calced_lines(i) = lineData(i,2)*lineData(i,4)*level_table(lineData(i,3), Te_best_i, ne_i, n0_best_i);
                end
                    
                sums_calc_best = [
                    calced_lines(1);
                    calced_lines(2);
                    calced_lines(3);
                    calced_lines(4);
                    calced_lines(5);
                    calced_lines(6) + calced_lines(7);
                    calced_lines(8);
                    calced_lines(9) + calced_lines(10);
                    calced_lines(11);
                    calced_lines(12) + calced_lines(13);
                    calced_lines(14) + calced_lines(15);
                    calced_lines(16);
                    calced_lines(17) + calced_lines(18);
                    calced_lines(19);
                    calced_lines(20);
                    calced_lines(21);
                    calced_lines(22);
                    calced_lines(23);
                    calced_lines(24);
                    calced_lines(25)];

