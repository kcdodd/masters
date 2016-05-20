function [ haszero,  zero_pts] = getZeroPts( f00, f10, f01, f11 )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

            haszero = false;
            zero_i = 0;
            zero_pts = zeros(2, 2);
            s00 = sign(f00);
            s01 = sign(f01);
            s10 = sign(f10);
            s11 = sign(f11);
            
            if (s00 ~= s01)
                zero_i = zero_i + 1;
                
                zero_pts(zero_i, 1) = 0;
                zero_pts(zero_i, 2) = 1/(abs(f01/f00) + 1);
            end
            
            if (s10 ~= s11)
                zero_i = zero_i + 1;
                
                zero_pts(zero_i, 1) = 1;
                zero_pts(zero_i, 2) = 1/(abs(f11/f10) + 1);
            end
            
            if (s00 ~= s10)
                zero_i = zero_i + 1;
                
                zero_pts(zero_i, 1) = 1/(abs(f10/f00) + 1);
                zero_pts(zero_i, 2) = 0;
            end
            
            if (s01 ~= s11)
                zero_i = zero_i + 1;
                
                zero_pts(zero_i, 1) = 1/(abs(f11/f01) + 1);
                zero_pts(zero_i, 2) = 1;
            end
            
            if (zero_i > 2)
                disp('too many zeros');
            end
            
            if (zero_i == 1)
                disp('too few zeros');
            end
            
            if (zero_i == 2)
                haszero = true;
            end
end

