function [sol] = loadRadialIonSolution(filename)
    
    fid = fopen(filename);
    
        tmp = textscan(fid, '%f', 1);
        nr = tmp{1,1};
        
        tmp = textscan(fid, '%f', 1);
        nz = tmp{1,1};
        
        tmp = textscan(fid, '%f', 1);
        nrbd = tmp{1,1};
        
        tmp = textscan(fid, '%f', 1);
        nzbd = tmp{1,1};
        
        tmp = textscan(fid, '%f', 1);
        nv = tmp{1,1};
        
        tmp = textscan(fid, '%f', 1);
        avg_n = tmp{1,1};
        
        tmp = textscan(fid, '%f', 1);
        avg_ion_rate = tmp{1,1};
    
    n = zeros(nr, nz);
    vr = zeros(nr, nz);
    vz = zeros(nr, nz);
    
    nbd_top = zeros(nrbd,1);
    nbd_bottom = zeros(nrbd,1);
    nbd_inner = zeros(nzbd,1);
    nbd_outer = zeros(nzbd,1);
        
        for iz = 1:nzbd
                tmp = textscan(fid, '%f', 1);
                nbd_outer(iz) = tmp{1,1};
        end

        for iz = 1:nzbd
                tmp = textscan(fid, '%f', 1);
                nbd_inner(iz) = tmp{1,1};
        end
        
        for ir = 1:nrbd
                tmp = textscan(fid, '%f', 1);
                nbd_top(ir) = tmp{1,1};
        end
        
        for ir = 1:nrbd
                tmp = textscan(fid, '%f', 1);
                nbd_bottom(ir) = tmp{1,1};
        end
        
        for iz = 1:nz
            for ir = 1:nr
                tmp = textscan(fid, '%f', 1);
                n(ir, iz) = tmp{1,1};
            end
        end

    

        for iz = 1:nz
            for ir = 1:nr
                tmp = textscan(fid, '%f', 1);
                vr(ir, iz) = tmp{1,1};
            end
        end

    

        for iz = 1:nz
            for ir = 1:nr
                tmp = textscan(fid, '%f', 1);
                vz(ir, iz) = tmp{1,1};
            end
        end

    
    fclose(fid);
    
    n_sum = sum(n(2:(nr-1),1:(nz-1)),2)/(nz-1);
    
    
    sol = struct('n', n, 'vr', vr, 'vz', vz, 'n_sum', n_sum, 'avg_n', avg_n, 'avg_ion_rate', avg_ion_rate);
end