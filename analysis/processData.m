function [processedData]=processData(rawdata, calibration, bg, peakratio, lineintegral, culldata, linedata, leveldata, iongauge_factor)

    num = size(rawdata,1);

    num_lines = size(linedata, 1);
    num_levels = size(leveldata,1);
    
    radiance = cell(num,1);
    npop = cell(num,1);
    nArI = cell(num,1);
    nArII = cell(num,1);
    
    L_lowpass = cell(num,1);
    L_highpass = cell(num,1);
    
    P_avg = cell(num,1);
    n0 = cell(num,1);
    integrations = cell(num,1);
    r = cell(num,1);
    
    spec_radiance_cells = cell(num,1);
    spec_counts_cells = cell(num,1);
   
    T0 = 300;
    kT0 = 1.38e-23*T0;

    for i = 1:num
        %average spec data
        numscans = size(rawdata(i,1).spec,2);
        numwl = size(rawdata(i,1).spec,1);

        % figure out how many different int times there are

        ints = unique(rawdata(i,1).int);
        integrations{i,1}=ints;
        num_ints = size(ints,1);

        spec_avg = zeros(numwl, num_ints);
        num_avged = zeros(num_ints);

        max_t = rawdata(i,1).t_spec(numscans) - culldata.t_spec_end;

        % average desired spec
        for j = 1:numscans
            if (rawdata(i,1).t_spec(j) >= culldata.t_spec_start && rawdata(i,1).t_spec(j) <= max_t)
                el = find(ints == rawdata(i,1).int(j));
                % subtract background when adding to the average, but don't
                % make negative
                spec_avg(:,el) = spec_avg(:,el) + max(rawdata(i,1).spec(:,j) - bg, 0);
                num_avged(el) = num_avged(el) + 1;
            end
        end

        % divide by number of scans used for each integration time
        for j = 1:num_ints
            
            spec_avg(:,j) = spec_avg(:,j)/num_avged(j);
            
        end
        
        spec_counts_cells = spec_avg;

        spec_radiance = zeros(size(spec_avg));

        % calibrate spec data
        for j = 1:num_ints
            spec_radiance(:,j) = spec_avg(:,j)./(1e-3*ints(j)*calibration);
        end

        spec_radiance_cells{i,1} = spec_radiance;
        
        % integrate lines
        
        radiance{i, 1} = zeros(num_lines, num_ints);
        npop{i, 1} = zeros(num_levels,1);
        npop_N = zeros(num_levels,1);

        % search for each line
        for j = 1:num_lines

            % each integration time
            for k = 1:num_ints
                peak = 0;
                for w = 1:numwl
                    if (rawdata(i,1).wl(w) > linedata(j, 1) - linedata(j, 2) && rawdata(i,1).wl(w) < linedata(j, 1) + linedata(j, 2))
                        if (spec_radiance(w, k) > peak)
                            peak = spec_radiance(w, k);
                        end
                    end
                end
                
                radiance{i, 1}(j, k) = peakratio*peak;
                

                
                for l = 1:num_levels
                    if (linedata(j,3) == leveldata(l,1) && linedata(j,4) == leveldata(l,2) && linedata(j,6) == ints(k))
                        npop{i, 1}(l) = npop{i, 1}(l) + radiance{i, 1}(j, k)/(lineintegral*linedata(j,5));
                        npop_N(l) = npop_N(l) + 1;
                        break;
                    end
                end
                
            end
        end
        
        npop{i, 1} = npop{i, 1}./npop_N;
        
        nArI{i,1} = 34*(npop{i, 1}(1)/3 + npop{i, 1}(2)/20 + npop{i, 1}(3)/8 + npop{i, 1}(4)/3)/4;
        nArII{i,1} = 42*(npop{i, 1}(5)/12 + npop{i, 1}(6)/20 + npop{i, 1}(7)/10)/3;
        
        % determin low/high pass sums
        %low-pass range
        lowpass_wl_min = 415;
        lowpass_wl_max = 515;
        
        %high-pass range
        highpass_wl_min = 575;
        highpass_wl_max = 1600;
        
        L_lowpass{i,1} = zeros(num_ints, 1);
        L_highpass{i,1} = zeros(num_ints, 1);
        
        for k = 1:num_ints
            lowsum = 0;
            highsum = 0;
            
            for w = 1:numwl
                
                if (rawdata(i,1).wl(w) > lowpass_wl_min && rawdata(i,1).wl(w) < lowpass_wl_max)
                    lowsum = lowsum + spec_radiance(w, k);
                end

                if (rawdata(i,1).wl(w) > highpass_wl_min && rawdata(i,1).wl(w) < highpass_wl_max)
                    highsum = highsum + spec_radiance(w, k);
                end
            end
            
            L_lowpass{i,1}(k) = lowsum;
            L_highpass{i,1}(k) = highsum;
        end
        
        %average desired pressure

        numscans = size(rawdata(i,1).P,1);

        P_avg{i,1} = 0;
        num_avged = 0;

        for j = 1:numscans
            if (rawdata(i,1).t_P(j) <= culldata.t_P_start)
                P_avg{i,1} = P_avg{i,1} + rawdata(i,1).P(j);
                num_avged = num_avged+1;
            end
        end

        P_avg{i,1} = P_avg{i,1}/num_avged;
        
        P_avg{i,1} = P_avg{i,1}/iongauge_factor;
        
        n0{i,1} = 133.3*P_avg{i,1}/kT0;
        
        r{i,1} = rawdata(i,1).r;
    end

    
    processedData = struct('nArI', nArI, 'nArII', nArII, 'npop', npop, 'leveldata', leveldata, 'linedata', linedata, 'spec_counts', spec_counts_cells, 'spec_radiance', spec_radiance_cells, 'radiance', radiance, 'L_lowpass', L_lowpass, 'L_highpass', L_highpass, 'pressure', P_avg, 'n0', n0, 'int', integrations, 'r', r);
end