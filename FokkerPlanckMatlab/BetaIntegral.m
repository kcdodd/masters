function [ I, b] = BetaIntegral(beta_max, N)
%BINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

    I=zeros(N,1);
    b=zeros(N,1);
    
    dbeta = beta_max/(N-1);
    
    beta = 0;

    for i=1:N
        b(i) = beta;
        
        f=@(t) (beta*sin(t)./sqrt(1+beta^2-2*beta*cos(t))).^3;
        
        I(i) = quadgk(f, 0, pi);
        
        beta = beta + dbeta;
    end
end

