function [ corrections ] = load_evdf_solutions( dirname, n0, num_velocities)
%LOAD_EVDF_SOLUTIONS Summary of this function goes here
%   Detailed explanation goes here

    constants;
    
    Te = [0.5, 1, 2, 5, 10, 20];
    ne = [1E16, 2E16, 5E16, 1E17, 2E17, 5E17, 1E18, 2E18];
    

    num_Te = size(Te, 2);
    num_ne = size(ne, 2);
    num_n0 = size(n0, 2);

    Ar_I_exc_corrections = ones(num_Te, num_ne, num_n0);
    Ar_I_ion_corrections = ones(num_Te, num_ne, num_n0);
    Ar_II_exc_corrections = ones(num_Te, num_ne, num_n0);
    
    files = dir(dirname);
    num_files=size(files,1);
    
        for i = 1:num_files

            if ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')
                    filename = [dirname  '\'  files(i).name];

                    fid = fopen(filename);

                    tmp = textscan(fid,'%f',1);
                    Te_i = tmp{1,1};
                    tmp = textscan(fid,'%f',1);
                    ne_i = tmp{1,1};
                    tmp = textscan(fid,'%f',1);
                    n0_i = tmp{1,1};
                    
                    tmp = textscan(fid,'%f',1);
                    tmp = textscan(fid,'%f',1);
                    
                    tmp = textscan(fid,'%f',1);
                    dv_i = tmp{1,1};
                    
                    f = zeros(num_velocities,1);
                    
                    for j = 1:num_velocities
                        tmp = textscan(fid,'%f',1);
                        f(j) = tmp{1,1};
                    end
                    
                    fclose(fid);
                    
                    Te_index = 0;
                    ne_index = 0;
                    n0_index = 0;
                    
                    for j = 1:num_Te
                        
                        if (Te(j) == Te_i)
                            Te_index = j;
                            break;
                        end
                    end
                    
                    for j = 1:num_ne
                        
                        if (ne(j) == ne_i)
                            ne_index = j;
                            break;
                        end
                    end
                    
                    for j = 1:num_n0
                        
                        if (n0(j) == n0_i)
                            n0_index = j;
                            break;
                        end
                    end
                    
                    if (Te_index == 0)
                        % not found
                        disp(['index not found Te=', num2str(Te_i)]);
                        break;
                    end
                    
                    if (ne_index == 0)
                        % not found
                        disp(['index not found ne=', num2str(ne_i)]);
                        break;
                    end
                    
                    if (n0_index == 0)
                        % not found
                        disp(['index not found n0=', num2str(n0_i)]);
                        break;
                    end
                    
                    
                    v = ((0:(num_velocities-1))')*dv_i;
                    
                  
                    
                    f_0 = exp(-v.^2*0.5*const_me/(const_e*Te_i));
                    f_0 = f_0/calc_integral_total((4*pi*v.^2).*f_0, dv_i);
                    


                    Ar_I_exc_corrections(Te_index, ne_index, n0_index) = min(1, calc_correction_factor(11.5, f, f_0, v));
                    Ar_I_ion_corrections(Te_index, ne_index, n0_index) = min(1, calc_correction_factor(15.76, f, f_0, v));
                    Ar_II_exc_corrections(Te_index, ne_index, n0_index) = min(1, calc_correction_factor(19.2, f, f_0, v));

                    
            end
            
        end
        %{
        for n = 1: num_n0
            
            for r = 1:num_Te
                for c = 1:num_ne

                    if (r < 3 || r > 5 || c > 5)
                        
                        r_i = r;
                        c_i = c;
                        
                        if ( r < 3)
                            r_i = 3;
                        end
                        
                        if ( r > 5)
                            r_i = 5;
                        end
                        
                        if (c > 5)
                            c_i = 5;
                        end
                        
                        Ar_I_exc_corrections(r, c, n) = Ar_I_exc_corrections(r_i, c_i, n);
                        Ar_I_ion_corrections(r, c, n) = Ar_I_ion_corrections(r_i, c_i, n);
                        Ar_II_exc_corrections(r, c, n) = Ar_II_exc_corrections(r_i, c_i, n);
                    end
                end

            end
        end
        %}
        
        corrections = struct('Te', Te, 'ne', ne, 'n0', n0, 'Ar_I_exc_corrections', Ar_I_exc_corrections, 'Ar_I_ion_corrections', Ar_I_ion_corrections, 'Ar_II_exc_corrections', Ar_II_exc_corrections);
end

