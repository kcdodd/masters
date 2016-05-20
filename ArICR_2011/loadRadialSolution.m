function [n, vr, vz, n0, vr0, vz0] = loadRadialSolution(filename, np, nr, nz)

    n = zeros(np, nr, nz);
    vr = zeros(np, nr, nz);
    vz = zeros(np, nr, nz);
    
    n0 = zeros(nr, nz);
    vr0 = zeros(nr, nz);
    vz0 = zeros(nr, nz);
    
    fid = fopen(filename);
    
    for ip = 1:np
        for iz = 1:nz
            for ir = 1:nr
                tmp = textscan(fid, '%f', 1);
                n(ip, ir, iz) = tmp{1,1};
            end
        end
    end
    
    for ip = 1:np
        for iz = 1:nz
            for ir = 1:nr
                tmp = textscan(fid, '%f', 1);
                vr(ip, ir, iz) = tmp{1,1};
            end
        end
    end
    
    for ip = 1:np
        for iz = 1:nz
            for ir = 1:nr
                tmp = textscan(fid, '%f', 1);
                vz(ip, ir, iz) = tmp{1,1};
            end
        end
    end
    
    fclose(fid);
    
    for ip = 1:np
        for iz = 1:nz
            for ir = 1:nr
                n0(ir, iz) = n0(ir, iz) + n(ip, ir, iz);
                vr0(ir, iz) = vr0(ir, iz) + n(ip, ir, iz)*vr(ip, ir, iz);
                vz0(ir, iz) = vz0(ir, iz) + n(ip, ir, iz)*vz(ip, ir, iz);
            end
        end
    end
    
        for iz = 1:nz
            for ir = 1:nr
                vr0(ir, iz) = vr0(ir, iz)/n0(ir, iz);
                vz0(ir, iz) = vz0(ir, iz)/n0(ir, iz);
            end
        end
    
end