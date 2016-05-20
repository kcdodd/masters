function [ derivative ] = calc_derivative(y_i, dx)
%CALC_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here

    N = size(y_i, 1);

    derivative = zeros(N,1);
    
    for i = 3:(N-2)
        derivative(i) = (-y_i(i+2) + 8*y_i(i+1) - 8*y_i(i-1) + y_i(i-2))/(12*dx);
    end
    
    derivative(1) = 0;
    derivative(2) = (-y_i(4) + 8*y_i(3) - 8*y_i(1) + y_i(2))/(12*dx);
    
    
    derivative(N-1) = (3*y_i(N) + 10*y_i(N-1) -18*y_i(N-2) + 6*y_i(N-3) - y_i(N-4))/(12*dx);
    derivative(N) = (25*y_i(N) - 48*y_i(N-1) + 36*y_i(N-2) - 16*y_i(N-3) + 3*y_i(N-4))/(12*dx);
    
    %derivative = smooth(derivative);
end

