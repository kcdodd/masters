function [Te, dTe, ne, dne, r, n0_avg, n0_ratio] = calcNeTe_3(level_table_ArI, level_table_ArII, processedData, ArI_error, ArII_error, ne_pkratio, vol_ratio, inelastic_corrections, n0_min_ratio)
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
    
    n0_interp_i = 0;
    n0_interp_factor = 0;
    
    for i = 1:size(inelastic_corrections.n0,2)
        
        if (inelastic_corrections.n0(i) >= n0_avg)
            n0_interp_i = i-1;
            n0_interp_factor = (n0_avg - inelastic_corrections.n0(i-1))/(inelastic_corrections.n0(i)-inelastic_corrections.n0(i-1));
            break;
        end
    
    end
    
    
    
    % sum over all levels to create effective total population of
    % measurable excited levels
    
    if fulltables
    
        tableArI_tot = 34*(level_table_ArI(1,1).npop/3 + level_table_ArI(2,1).npop/20 + level_table_ArI(3,1).npop/8 + level_table_ArI(4,1).npop/3)/4;

        tableArII_tot = 42*(level_table_ArII(1,1).npop/12 + level_table_ArII(2,1).npop/20 + level_table_ArII(3,1).npop/10)/3;
        
    else
        tableArI_tot = 31*(level_table_ArI(1,1).npop/3 + level_table_ArI(2,1).npop/20 + level_table_ArI(3,1).npop/8)/3;

        tableArII_tot = 32*(level_table_ArII(1,1).npop/12 + level_table_ArII(2,1).npop/20)/2;
    end
    
    % correct level table populations for non-maxwellian distribution from inelastic
    % collisions.
   
    tableArI_tot = tableArI_tot.*((1-n0_interp_factor)*inelastic_corrections.Ar_I_exc_corrections(:,:,n0_interp_i) + n0_interp_factor*inelastic_corrections.Ar_I_exc_corrections(:,:,n0_interp_i+1));
    tableArII_tot = tableArII_tot.*((1-n0_interp_factor)*inelastic_corrections.Ar_II_exc_corrections(:,:,n0_interp_i) + n0_interp_factor*inelastic_corrections.Ar_II_exc_corrections(:,:,n0_interp_i+1));

    
    for u = 1:100

        ni_avg = 0;
        
        n0_ratio_old = n0_ratio;

        for i = 1:n_r
            r(i) = processedData(i,1).r;


            nes = interp1(log10(level_table_ArI(1,1).ne),1:0.0625:size(level_table_ArI(1,1).ne,2));
            Tes = interp1(level_table_ArI(1,1).Te, 1:0.0625:size(level_table_ArI(1,1).Te,2));
            
            if fulltables
                npopArI_tot = 34*(processedData(i,1).npop(1)/3 + processedData(i,1).npop(2)/20 + processedData(i,1).npop(3)/8 + processedData(i,1).npop(4)/3)/(4*n0_min_ratio(i)*n0_ratio*processedData(i,1).n0);
                
                npopArII_tot = 42*(processedData(i,1).npop(5)/12 + processedData(i,1).npop(6)/20 + processedData(i,1).npop(7)/10)/3;
            else
                npopArI_tot = 31*(processedData(i,1).npop(1)/3 + processedData(i,1).npop(2)/20 + processedData(i,1).npop(3)/8)/(3*n0_min_ratio(i)*n0_ratio*processedData(i,1).n0);
                
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
        
        n0_ratio = (n0_avg - vol_ratio*ni_avg)/n0_avg;
        
        %n0_ratio = 0.9;
        
        if (abs((n0_ratio - n0_ratio_old)*2/(n0_ratio + n0_ratio_old)) < err)
            % stable n0_ratio found
            break;
        end
        
    end
    
    %calculate uncertainties
        for i = 1:n_r

            nes = interp1(log10(level_table_ArI(1,1).ne),1:0.0625:size(level_table_ArI(1,1).ne,2));
            Tes = interp1(level_table_ArI(1,1).Te, 1:0.0625:size(level_table_ArI(1,1).Te,2));
            
            if fulltables
                npopArI_tot = 34*(processedData(i,1).npop(1)/3 + processedData(i,1).npop(2)/20 + processedData(i,1).npop(3)/8 + processedData(i,1).npop(4)/3)/(4*n0_min_ratio(i)*n0_ratio*processedData(i,1).n0);
                
                npopArII_tot = 42*(processedData(i,1).npop(5)/12 + processedData(i,1).npop(6)/20 + processedData(i,1).npop(7)/10)/3;
            else
                npopArI_tot = 31*(processedData(i,1).npop(1)/3 + processedData(i,1).npop(2)/20 + processedData(i,1).npop(3)/8)/(3*n0_min_ratio(i)*n0_ratio*processedData(i,1).n0);
                
                npopArII_tot = 32*(processedData(i,1).npop(5)/12 + processedData(i,1).npop(6)/20)/2;
            end
            
            fI = interp2(log10(tableArI_tot) - log10(npopArI_tot), 4,'cubic');
            fII = interp2(log10(tableArII_tot) - log10(npopArII_tot), 4,'cubic');
            
            npopArI_tot_high = npopArI_tot*(1 + ArI_error);
            npopArI_tot_low = npopArI_tot*(1 - ArI_error);

            npopArII_tot_high = npopArII_tot*(1 + ArII_error);
            npopArII_tot_low = npopArII_tot*(1 - ArII_error);
            
            fI_high = interp2(log10(tableArI_tot) - log10(npopArI_tot_high), 4,'cubic');
            fI_low = interp2(log10(tableArI_tot) - log10(npopArI_tot_low), 4,'cubic');
            
            fII_high = interp2(log10(tableArII_tot) - log10(npopArII_tot_high), 4,'cubic');
            fII_low = interp2(log10(tableArII_tot) - log10(npopArII_tot_low), 4,'cubic');

            [ne_1, Te_1] = findZeroIntersection(nes, Tes, fI_high, fII);
            [ne_2, Te_2] = findZeroIntersection(nes, Tes, fI_low, fII);
            
            dTe_I = abs(Te_1-Te_2)/2;
            dne_I = abs(10^ne_1-10^ne_2)/2;
            
            [ne_1, Te_1] = findZeroIntersection(nes, Tes, fI, fII_high);
            [ne_2, Te_2] = findZeroIntersection(nes, Tes, fI, fII_low);
            
            dTe_II = abs(Te_1-Te_2)/2;
            dne_II = abs(10^ne_1-10^ne_2)/2;
            
            dTe(i) = sqrt(dTe_I^2 + dTe_II^2);
            dne(i) = sqrt(dne_I^2 + dne_II^2);
        end


end

