function [ out ] = f1(x)
%F1 Summary of this function goes here
%   Detailed explanation goes here

    g = @(t) exp(x*(1-t))./t;
    
    out = quadgk(g, 1, 100);
end

