function [ out ] = F(x, A, B, C, D)
%F Summary of this function goes here
%   Detailed explanation goes here
A_i = A*(1-x*f1(x));
B_i = B*(1 + x - x*(2+x)*f1(x));
C_i = C*f1(x);
D_i = D*x*f2(x);

out = A_i + B_i + C_i + D_i;

end

