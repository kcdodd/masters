function [ r ] = calc_delta2(y1, y2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N1 = max(size(y1, 1), size(y1, 2));
N2 = max(size(y2, 1), size(y2, 2));

if (N1 ~= N2)
    disp('functions must be the same length');
end

y1_sum = 0;
y2_sum = 0;
y1y2_sum = 0;
y1_2_sum = 0;
y2_2_sum = 0;

for i = 1:N1
    y1_sum = y1_sum + y1(i);
    y2_sum = y2_sum + y2(i);
    y1y2_sum = y1y2_sum + y1(i)*y2(i);
    y1_2_sum = y1_2_sum + y1(i)^2;
    y2_2_sum = y2_2_sum + y2(i)^2;
    
end

avg_y1 = y1_sum/N1;
avg_y2 = y2_sum/N1;

s1 = sqrt((y1_2_sum - N1*avg_y1^2)/(N1-1));
s2 = sqrt((y2_2_sum - N1*avg_y2^2)/(N1-1));

r = (y1y2_sum - N1*avg_y1*avg_y2)/(s1*s2*(N1-1));


end

