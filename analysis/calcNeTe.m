function [Te, dTe, ne, dne,  r, n0_ratio] = calcNeTe(level_table_ArI, level_table_ArII, processedData)
%CALCNETE Summary of this function goes here
%   Detailed explanation goes here

    use_ArI_7 = true;

    n_r = size(processedData, 1);
    
    Te = zeros(n_r,1);
    dTe = zeros(n_r,1);
    ne = zeros(n_r,1);
    dne = zeros(n_r,1);
    r = zeros(n_r,1);
    
    n0_ratio = 1;
    n0_avg = 0;
    
    for i = 1:n_r
        n0_avg = n0_avg + processedData(i,1).n0;
    end
    
    n0_avg = n0_avg/n_r;
    
    for u = 1:5

        ni_avg = 0;

        for i = 1:n_r
            r(i) = processedData(i,1).r;


            nes = interp1(log10(level_table_ArI(1,1).ne),1:0.0625:size(level_table_ArI(1,1).ne,2));
            Tes = interp1(level_table_ArI(1,1).Te, 1:0.0625:size(level_table_ArI(1,1).Te,2));


            % A - 1
            nI_reduc = processedData(i,1).npop(2)/(n0_ratio*processedData(i,1).n0);
            fI_7 = interp2(log10(level_table_ArI(2,1).npop) - log10(nI_reduc), 4,'cubic');

            nII = processedData(i,1).npop(5);
            fII_13 = interp2(log10(level_table_ArII(1,1).npop) - log10(nII), 4,'cubic');

            [ne_A1, Te_A1] = findZeroIntersection(nes, Tes, fI_7, fII_13);

            % A - 2
            nII = processedData(i,1).npop(6);
            fII_14 = interp2(log10(level_table_ArII(2,1).npop) - log10(nII), 4,'cubic');

            [ne_A2, Te_A2] = findZeroIntersection(nes, Tes, fI_7, fII_14);


            % B - 1

            nI_reduc = processedData(i,1).npop(1)/(n0_ratio*processedData(i,1).n0);
            fI_6 = interp2(log10(level_table_ArI(1,1).npop) - log10(nI_reduc), 4,'cubic');

            [ne_B1, Te_B1] = findZeroIntersection(nes, Tes, fI_6, fII_13);

            % B - 2

            nI_reduc = processedData(i,1).npop(3)/(n0_ratio*processedData(i,1).n0);
            fI_8 = interp2(log10(level_table_ArI(3,1).npop) - log10(nI_reduc), 4,'cubic');

            [ne_B2, Te_B2] = findZeroIntersection(nes, Tes, fI_8, fII_13);

            % B - 3

            [ne_B3, Te_B3] = findZeroIntersection(nes, Tes, fI_6, fII_14);

            % B - 4

            [ne_B4, Te_B4] = findZeroIntersection(nes, Tes, fI_8, fII_14);

            if (use_ArI_7)
                
                %Te_max = max([Te_A1 Te_A2 Te_B1 Te_B2 Te_B3 Te_B4]);
                %Te_min = min([Te_A1 Te_A2 Te_B1 Te_B2 Te_B3 Te_B4]);

                %ne_max = max([ne_A1 ne_A2 ne_B1 ne_B2 ne_B3 ne_B4]);
                %ne_min = min([ne_A1 ne_A2 ne_B1 ne_B2 ne_B3 ne_B4]);

                Te_max = max([Te_A1 Te_A2]);
                Te_min = min([Te_A1 Te_A2]);

                ne_max = max([ne_A1 ne_A2]);
                ne_min = min([ne_A1 ne_A2]);
            else
                
                Te_max = max([Te_B1 Te_B2 Te_B3 Te_B4]);
                Te_min = min([Te_B1 Te_B2 Te_B3 Te_B4]);

                ne_max = max([ne_B1 ne_B2 ne_B3 ne_B4]);
                ne_min = min([ne_B1 ne_B2 ne_B3 ne_B4]);
            
            end

            %compute average and error

            ne(i) = (10^(ne_max) + 10^(ne_min))/2;
            dne(i) = abs(10^(ne_max) - 10^(ne_min))/2;
            Te(i) = (Te_max + Te_min)/2;
            dTe(i) = abs(Te_max - Te_min)/2;
            
            ni_avg = ni_avg + ne(i);
        end
        
        ni_avg = ni_avg/n_r;
        
        n0_ratio = (n0_avg - ni_avg)/n0_avg;
    end

end

