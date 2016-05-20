function Te = calc_Te (wavelengths, emission, ne)

Te_l = [0.5 1 2 5 10 20];
ne_l = [1e10 2e10 5e10 1e11 2e11 5e11 1e12 2e12];

if (ne < min(ne_l) || ne > max(ne_l))
	disp('density outside bounds');
	return
end

hc = 6.626e-34*3e8;

% center wavelength, level, factor, low, high
% level 1=13, 2=14, 3=15
lineData = [
	434 2 hc*(5.74e7*4/(433.12e-9)+1.171e8*8/(434.81e-9))/20 2 2;
	442.8 2 hc*(8.17e7*6/(442.6e-9)+5.69e7*4/(443.02e-9))/20 2 2;
	480.6 1 hc*(7.80e7*6/(480.6e-9))/12 2 2;
	487.98 3 hc*(8.23e7*6/(487.98e-9))/10 2 2];
	
% level 13 rate data
data_1 = [
	3.23e-23	7.89e-23	2.99e-22	8.98e-22	2.79e-21	1.21e-20	3.45e-20	9.48e-20;
	4.34e-15	1.03e-14	3.66e-14	1.04e-13	3.01e-13	1.18e-12	3.18e-12	8.37e-12;
	4.29e-11	9.77e-11	3.24e-10	8.65e-10	2.36e-09	8.58e-09	2.20e-08	5.62e-08;
	8.69e-09	1.91e-08	5.90e-08	1.48e-07	3.83e-07	1.32e-06	3.30e-06	8.33e-06;
	4.08e-08	8.79e-08	2.63e-07	6.42e-07	1.63e-06	5.54e-06	1.38e-05	3.49e-05;
	6.52e-08	1.39e-07	4.02e-07	9.77e-07	2.45e-06	8.35e-06	2.09e-05	5.31e-05];

%level 14 rate data
data_2 = [
	1.83e-23	4.08e-23	1.33e-22	3.70e-22	1.15e-21	5.65e-21	1.88e-20	5.89e-20;
	3.90e-15	8.46e-15	2.59e-14	6.68e-14	1.88e-13	8.05e-13	2.47e-12	7.42e-12;
	4.75e-11	1.01e-10	2.95e-10	7.15e-10	1.86e-09	7.12e-09	2.04e-08	5.88e-08;
	1.06e-08	2.22e-08	6.19e-08	1.43e-07	3.48e-07	1.22e-06	3.33e-06	9.32e-06;
	5.04e-08	1.04e-07	2.86e-07	6.44e-07	1.53e-06	5.23e-06	1.40e-05	3.89e-05;
	7.99e-08	1.65e-07	4.45e-07	9.90e-07	2.32e-06	7.79e-06	2.07e-05	5.77e-05];

%level 15 rate data
data_3 = [
	2.48e-23	5.01e-23	1.30e-22	2.75e-22	6.24e-22	2.15e-21	6.16e-21	1.85e-20;
	8.59e-15	1.72e-14	4.29e-14	8.65e-14	1.78e-13	5.06e-13	1.23e-12	3.24e-12;
	1.33e-10	2.64e-10	6.50e-10	1.27e-09	2.50e-09	6.32e-09	1.38e-08	3.24e-08;
	3.21e-08	6.37e-08	1.56e-07	3.02e-07	5.79e-07	1.38e-06	2.78e-06	6.02e-06;
	1.56e-07	3.09e-07	7.56e-07	1.47e-06	2.80e-06	6.55e-06	1.28e-05	2.67e-05;
	2.62e-07	5.21e-07	1.28e-06	2.49e-06	4.76e-06	1.11e-05	2.13e-05	4.32e-05];
	
sums_spec = zeros(size(lineData, 1), 1);
lines = zeros(size(lineData, 1), 1);
    

for i = 1:size(lineData, 1)
    lines(i) = lineData(i, 1);
    
    for w = 1:size(wavelengths, 1)
        if (wavelengths(w) > lineData(i, 1) - lineData(i, 4) && wavelengths(w) < lineData(i, 1) + lineData(i, 5))
            sums_spec(i) = sums_spec(i) + emission(w);
        end
    end
    
end

level_rate = zeros(3,1);

level_rate(1) = (sums_spec(3)/lineData(3,3))/ne;
level_rate(2) = ((sums_spec(1)/lineData(1,3) + sums_spec(2)/lineData(2,3))/2)/ne;
level_rate(3) = (sums_spec(4)/lineData(4,3))/ne;

disp(level_rate);

% calc ne interp factor
ne_factor = 0;
ne_i = 0;
for i = 1:5
	if (ne > ne_l(i) && ne <= ne_l(i+1))
		ne_i = i;
		ne_factor = (ne - ne_l(i))/(ne_l(i+1) - ne_l(i));
	end
end

% calc temp from 13
temp_1 = 0;

for i = 1:5
	rate_min = ne_factor*(data_1(i,ne_i+1) - data_1(i,ne_i)) + data_1(i,ne_i);
	rate_max = ne_factor*(data_1(i+1,ne_i+1) - data_1(i+1,ne_i)) + data_1(i+1,ne_i);
	
	if (level_rate(1) > rate_min && level_rate(1) < rate_max)
		temp_fact = (level_rate(1) - rate_min)/(rate_max - rate_min);
		temp_1 = temp_fact*(Te_l(i+1) - Te_l(i)) + Te_l(i);
	end
end

% calc temp from 14
temp_2 = 0;

for i = 1:5
	rate_min = ne_factor*(data_2(i,ne_i+1) - data_2(i,ne_i)) + data_2(i,ne_i);
	rate_max = ne_factor*(data_2(i+1,ne_i+1) - data_2(i+1,ne_i)) + data_2(i+1,ne_i);
	
	if (level_rate(2) > rate_min && level_rate(2) < rate_max)
		temp_fact = (level_rate(2) - rate_min)/(rate_max - rate_min);
		temp_2 = temp_fact*(Te_l(i+1) - Te_l(i)) + Te_l(i);
	end
end

% calc temp from 15
temp_3 = 0;

for i = 1:5
	rate_min = ne_factor*(data_3(i,ne_i+1) - data_3(i,ne_i)) + data_3(i,ne_i);
	rate_max = ne_factor*(data_3(i+1,ne_i+1) - data_3(i+1,ne_i)) + data_3(i+1,ne_i);
	
	if (level_rate(3) > rate_min && level_rate(3) < rate_max)
		temp_fact = (level_rate(3) - rate_min)/(rate_max - rate_min);
		temp_3 = temp_fact*(Te_l(i+1) - Te_l(i)) + Te_l(i);
	end
end

if (temp_1 == 0 || temp_2 == 0 || temp_3 == 0)
	disp(['Te could not be determined: ' num2str(temp_1) ', ' num2str(temp_2) ', ' num2str(temp_3) '.']);
	return
end

% average values
Te = (temp_1 + temp_2 + temp_3)/3;