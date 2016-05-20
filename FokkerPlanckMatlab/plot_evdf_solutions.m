function plot_evdf_solutions( dirname, Te, ne, n0, num_velocities )
%LOAD_EVDF_SOLUTIONS Summary of this function goes here
%   Detailed explanation goes here

    constants;

    num_Te = size(Te, 2);
    num_ne = size(ne, 2);
    num_n0 = size(n0, 2);

    files = dir(dirname);
    num_files=size(files,1);
    
    num_loaded = 0;
    
    hold on;
                    
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
                        %disp(['index not found Te=', num2str(Te_i)]);
                        continue;
                    end
                    
                    if (ne_index == 0)
                        % not found
                        %disp(['index not found ne=', num2str(ne_i)]);
                        continue;
                    end
                    
                    if (n0_index == 0)
                        % not found
                        %disp(['index not found n0=', num2str(n0_i)]);
                        continue;
                    end
                    
                    
                    v = ((0:(num_velocities-1))')*dv_i;
                    
                    %plot(v, (4*pi*v.^2).*f);
                    %plot(v, f);
                    plot(0.5*const_me*v.^2/const_e, log(f));
                    
                    num_loaded = num_loaded + 1;
                    
                    disp('plotted');
                     
            end
            
        end
        
                    f_0 = exp(-v.^2*0.5*const_me/(const_e*Te(1)));
                    f_0 = f_0/calc_integral_total((4*pi*v.^2).*f_0, v(2));
                    
                    %plot(v, (4*pi*v.^2).*f_0);
                    plot(0.5*const_me*v.^2/const_e, log(f_0));
                    %plot(v, f_0);
                    
                    hold off;
                    
                    %num_loaded
        
end

