function [ M, alphas ] = calc_M(max_alpha, num, outfile )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    
    
    alphas = (0:num)*max_alpha/num;
    M = zeros(size(alphas));
    
    for i = 0:num
        f = @(beta) (beta.^2).*exp(-beta.^2-alphas(i+1)./beta);
        M(i+1) = quadgk(f, 0, 1) + quadgk(f, 1, 5);
    end
    
    fileID = fopen(outfile, 'w');
    
    fprintf(fileID, '%u\r\n', (num+1));
    
    fprintf(fileID, '%f\r\n', max_alpha);
    
    for i = 0:num
        fprintf(fileID, '%E\r\n', M(i+1));
    end
    
    fclose(fileID);
end
