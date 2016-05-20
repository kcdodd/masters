function [lines_summed, line_table] = calcLineTable(ne, Te, n0, Nediv, xrad)

line_table = zeros(15, size(ne,2), size(Te,2), size(n0,2));

for ne_i = 1:size(ne,2)
    
    for Te_i = 1:size(Te,2)
        
        for n0_i = 1:size(n0,2)
            
            [xn0, pop, gain, loss] = arcrmodel_new3(ne(ne_i), Te(Te_i), n0(n0_i), Nediv, xrad);
            
            [lines_all, emission] = calcLineEmission(pop);
            
            [sums_cr, lines_summed] = calc_summed_emission_cr(emission);
            
            line_table(:, ne_i, Te_i, n0_i) = sums_cr;
        end
        
    end
    
end
