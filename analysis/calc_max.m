function [ max_x, a, b, c ] = calc_max( x, y )
%UNTITLED Summary of this function goes here
%   for three points, uses quadratic interpolation to find maximum of
%   function

k = (x(3)-x(2))/(x(2)-x(1));

a = (y(3) - y(2) - k*(y(2)-y(1)))/(x(3)^2 - x(2)^2 - k*(x(2)^2-x(1)^2));

b = (y(2) - y(1) - a*(x(2)^2-x(1)^2))/(x(2)-x(1));

c = y(1) - a*x(1)^2 - b*x(1);

max_x = -b/(2*a);

end

