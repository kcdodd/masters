function [ out ] = f2( x )
%F2 Summary of this function goes here
%   Detailed explanation goes here

    g = @(t) log(t).*exp(x*(1-t))./t;
    
    out = quadgk(g, 1, 100);
end

