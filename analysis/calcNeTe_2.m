function [Te, ne,  r, n0_avg, n0_ratio] = calcNeTe_2(level_table_ArI, level_table_ArII, processedData, ne_pkratio)
%CALCNETE Summary of this function goes here
%   Detailed explanation goes here

    fulltables = true;

    n_r = size(processedData, 1);
    
    err = 1E-6;
    
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
    
    if fulltables
    
        tableArI_tot = 34*(level_table_ArI(1,1).npop/3 + level_table_ArI(2,1).npop/20 + level_table_ArI(3,1).npop/8 + level_table_ArI(4,1).npop/3)/4;

        tableArII_tot = 42*(level_table_ArII(1,1).npop/12 + level_table_ArII(2,1).npop/20 + level_table_ArII(3,1).npop/10)/3;
        
    else
        tableArI_tot = 31*(level_table_ArI(1,1).npop/3 + level_table_ArI(2,1).npop/20 + level_table_ArI(3,1).npop/8)/3;

        tableArII_tot = 32*(level_table_ArII(1,1).npop/12 + level_table_ArII(2,1).npop/20)/2;
    end
    
    for u = 1:100

        ni_avg = 0;
        
        n0_ratio_old = n0_ratio;

        for i = 1:n_r
            r(i) = processedData(i,1).r;


            nes = interp1(log10(level_table_ArI(1,1).ne),1:0.0625:size(level_table_ArI(1,1).ne,2));
            Tes = interp1(level_table_ArI(1,1).Te, 1:0.0625:size(level_table_ArI(1,1).Te,2));
            
            if fulltables
                npopArI_tot = 34*(processedData(i,1).npop(1)/3 + processedData(i,1).npop(2)/20 + processedData(i,1).npop(3)/8 + processedData(i,1).npop(4)/3)/(4*n0_ratio*processedData(i,1).n0);
                
                npopArII_tot = 42*(processedData(i,1).npop(5)/12 + processedData(i,1).npop(6)/20 + processedData(i,1).npop(7)/10)/3;
            else
                npopArI_tot = 31*(processedData(i,1).npop(1)/3 + processedData(i,1).npop(2)/20 + processedData(i,1).npop(3)/8)/(3*n0_ratio*processedData(i,1).n0);
                
                npopArII_tot = 32*(processedData(i,1).npop(5)/12 + processedData(i,1).npop(6)/20)/2;
            end
            
            
            
            
            fI = interp2(log10(tableArI_tot) - log10(npopArI_tot), 4,'cubic');
            
            

            fII = interp2(log10(tableArII_tot) - log10(npopArII_tot), 4,'cubic');

            [ne_i, Te_i] = findZeroIntersection(nes, Tes, fI, fII);

            ne(i) = 10^(ne_i);
            Te(i) = Te_i;
            
            ni_avg = ni_avg + ne(i);
        end
        
        ni_avg = max(ne)/ne_pkratio;
        
        n0_ratio = (n0_avg - ni_avg)/n0_avg;
        
        %n0_ratio = 0.9;
        
        if (abs((n0_ratio - n0_ratio_old)*2/(n0_ratio + n0_ratio_old)) < err)
            % stable n0_ratio found
            break;
        end
        
    end

end

