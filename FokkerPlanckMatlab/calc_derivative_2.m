function [ derivative_2 ] = calc_derivative_2(y_i, dx)
%CALC_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here

    N = size(y_i, 1);

    derivative_2 = zeros(N,1);
    
    for i = 3:(N-2)
        derivative_2(i) = (-y_i(i+2) + 16*y_i(i+1) - 30*y_i(i) + 16*y_i(i-1) - y_i(i-2))/(12*dx*dx);
    end
    
    derivative_2(2) = (-y_i(4) + 16*y_i(3) - 30*y_i(2) + 16*y_i(1) - y_i(2))/(12*dx*dx);
    
    derivative_2(1) = (-y_i(3) + 16*y_i(2) - 30*y_i(1) + 16*y_i(2) - y_i(3))/(12*dx*dx);
    
    derivative_2(N-1) = (11*y_i(N) - 20*y_i(N-1) + 6*y_i(N-2) + 4*y_i(N-3) - y_i(N-4))/(12*dx*dx);
    
    derivative_2(N) = (35*y_i(N) - 104*y_i(N-1) + 114*y_i(N-2) - 56*y_i(N-3) + 11*y_i(N-4))/(12*dx*dx);
    
    
    %derivative_2 = smooth(derivative_2);
end

