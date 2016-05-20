function [ smoothed_y ] = smooth( y_i )
%SMOOTH Summary of this function goes here
%   Detailed explanation goes here

    nf = max(size(y_i, 1),size(y_i, 2));

    kernal_x = -2:0.5:2;
    kernal = exp(-kernal_x.^2);
    kernal = kernal/sum(kernal);
    
    nk = size(kernal,2);
    nk_m = (nk-1)/2;
    
    smoothed_y = zeros(size(y_i));
    
    for i = 1:nf
        
        for j = 1:nk
            
            ij = i-nk_m+j-1;
            
            if (ij < 1)
                ij = abs(ij)+1;
            end
            
            if (ij > nf)
                ij = 2*nf-ij;
            end
            
            smoothed_y(i) = smoothed_y(i) + y_i(ij)*kernal(j);
        end
        
    end
    
end

