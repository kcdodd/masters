function [ level_table_ArII ] = correct_ArII( level_table_ArII, n0, ion_table )
%CORRECT_ARII Summary of this function goes here
%   Detailed explanation goes here
        
        for m = 1:size(level_table_ArII(1,1).Te,2)
            A = ion_table.C2(m)/ion_table.C1(m);

            for n = 1:size(level_table_ArII(1,1).ne,2)
                
                f = @(u) u - A*level_table_ArII(1,1).ne(n)/(n0*(2-exp(u)));
                
                x = fzero(f, 1);
                
                level_table_ArII(1,1).npop(m, n) = level_table_ArII(1,1).npop(m, n)/(2*exp(x)-1);
                level_table_ArII(2,1).npop(m, n) = level_table_ArII(2,1).npop(m, n)/(2*exp(x)-1);
                level_table_ArII(3,1).npop(m, n) = level_table_ArII(3,1).npop(m, n)/(2*exp(x)-1);
            end
        end
end

