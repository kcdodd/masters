function [ r ] = calc_r( y1, y2 )
%Calculates r statistic
%   y1 = actual data
%   y2 = fitting function

N1 = max(size(y1, 1), size(y1, 2));
N2 = max(size(y2, 1), size(y2, 2));

if (N1 ~= N2)
    disp('functions must be the same length');
end

error_sum = 0;
avg_sum = 0;
y1_2_sum = 0;

for i = 1:N1
    error_sum = error_sum + (y1(i)-y2(i))^2;
    avg_sum = avg_sum + y1(i);
    y1_2_sum = y1_2_sum + y1(i)^2;
end

avg = avg_sum/N1;
r = 1 - error_sum/(y1_2_sum - N1*avg^2);


end

