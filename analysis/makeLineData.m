function [ linedata] = makeLineData( )
%MAKELINEDATA Summary of this function goes here
%   Detailed explanation goes here

hc = 6.626e-34*3e8;

% wavelenght, line width, table, level, radiance coeff., integration time
% to use (ms)
linedata = [
    434, 3, 2, 2, (hc/(4*pi))*(5.74e7*4/(433.12e-9)+1.171e8*8/(434.81e-9))/20, 1000;
	442.8, 3, 2, 2, (hc/(4*pi))*(8.17e7*6/(442.6e-9)+5.69e7*4/(443.02e-9))/20, 1000;
	480.6, 3, 2, 1, (hc/(4*pi))*(7.80e7*6/(480.6e-9))/12, 1000;
	487.98, 3, 2, 3, (hc/(4*pi))*(8.23e7*6/(487.98e-9))/10, 1000;
    696.5431, 3, 1, 4, (hc/(4*pi))*(6.39e+06/(696.5431e-9)), 1000;
    706.7218, 3, 1, 3, (hc/(4*pi))*(3.80e+06*5/(706.7218e-9))/8, 1000;
    738.3980, 3, 1, 3, (hc/(4*pi))*(8.47e+06*5/(738.3980e-9))/8, 1000;
    763.5106, 3, 1, 2, (hc/(4*pi))*(2.45e+07*5/(763.5106e-9))/20, 200;
    794.8176, 3, 1, 3, (hc/(4*pi))*(1.86e+07*3/(794.8176e-9))/8, 1000;
    801, 3, 1, 2, (hc/(4*pi))*(4.90e+06*5/(800.6157e-9) + 9.28e+06*5/(801.4786e-9))/20, 1000;
    811, 4, 1, 2, (hc/(4*pi))*(2.5e+07*3/(810.3693e-9) + 3.31e+07*7/(811.5311e-9))/20, 200;
    826.4522, 3, 1, 4, (hc/(4*pi))*(1.53e+07/(826.4522e-9)), 1000;
    852.1442, 3, 1, 3, (hc/(4*pi))*(1.39e+07*3/(852.1442e-9))/8, 1000;
    912.2967, 3, 1, 1, (hc/(4*pi))*(1.89e+7/(912.2967e-9)), 200;
    922.4499, 3, 1, 2, (hc/(4*pi))*(5.03e+6*5/(922.4499e-9))/20, 1000;
    965.7786, 3, 1, 1, (hc/(4*pi))*(5.43e6/(965.7786e-9)), 1000];

end

