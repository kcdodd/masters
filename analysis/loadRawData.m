function [rawdata]=loadRawData(dirname, shotlist)

    numshots = size(shotlist, 2);

    wl=cell(numshots,1);
    I=cell(numshots,1);
    int=cell(numshots,1);
    P=cell(numshots,1);
    r=cell(numshots,1);
    t_I=cell(numshots,1);
    t_P=cell(numshots,1);

    wl_node='\machine::top.usb650:calib:wve_spec';
    spec_node='\machine::top.usb650:raw_spectra';
    r_node='\machine::top.usb650:setup:su_maj_r1';
    t_I_node='\machine::top.usb650:framing:frm_start';
    int_node='\machine::top.usb650:framing:frm_int';

    pres_node='pressure:p2_ion_gauge:pressure';
    t_P_node='pressure:p2_ion_gauge:time';

    files = dir(dirname);
    N=size(files,1);

    cur_shot = 0;

    for s = shotlist
    cur_shot = cur_shot + 1;

        for i = 1:N

            if ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')

                C=textscan(files(i).name,'%u_%s');
                shot=C{1};
                %type=C{2}{1}(1:4);

                if shot == s

                    filename = [dirname  '\'  files(i).name];

                    fid = fopen(filename);

                    line = fgets(fid);

                    while ischar (line)

                        if ~isempty(strfind(line, wl_node))

                            datatype = textscan(fid,'%u',1);
                            numdims = textscan(fid,'%u',1);
                            numels = textscan(fid,'%u',1);

                            wl{cur_shot,1}=zeros(numels{1,1},1);

                            for j=1:numels{1,1}
                                el = textscan(fid,'%u',1);
                                wl{cur_shot,1}(j)= el{1,1};
                            end
                        elseif ~isempty(strfind(line, spec_node))
                            datatype = textscan(fid,'%u',1);
                            numdims = textscan(fid,'%u',1);
                            numels1 = textscan(fid,'%u',1);
                            numels2 = textscan(fid,'%u',1);

                            I{cur_shot,1}=zeros(numels1{1,1},numels2{1,1});

                            for j=1:numels1{1,1}
                                for k=1:numels2{1,1}
                                    el = textscan(fid,'%f',1);
                                    I{cur_shot,1}(j,k)= el{1,1};
                                end
                            end
                        elseif ~isempty(strfind(line, r_node))
                            datatype = textscan(fid,'%u',1);
                            numdims = textscan(fid,'%u',1);
                            numels = textscan(fid,'%u',1);

                            el = textscan(fid,'%f',1);

                            r{cur_shot,1} = el{1,1};
                        elseif ~isempty(strfind(line, t_I_node))
                            datatype = textscan(fid,'%u',1);
                            numdims = textscan(fid,'%u',1);
                            numels = textscan(fid,'%u',1);

                            t_I{cur_shot,1}=zeros(numels{1,1},1);

                            for j=1:numels{1,1}
                                el = textscan(fid,'%f',1);
                                t_I{cur_shot,1}(j)= el{1,1};
                            end
                        elseif ~isempty(strfind(line, int_node))
                            datatype = textscan(fid,'%u',1);
                            numdims = textscan(fid,'%u',1);
                            numels = textscan(fid,'%u',1);

                            int{cur_shot,1}=zeros(numels{1,1},1);

                            for j=1:numels{1,1}
                                el = textscan(fid,'%f',1);
                                int{cur_shot,1}(j)= el{1,1};
                            end
                        elseif ~isempty(strfind(line, pres_node))
                            datatype = textscan(fid,'%u',1);
                            numdims = textscan(fid,'%u',1);
                            numels = textscan(fid,'%u',1);

                            P{cur_shot,1}=zeros(numels{1,1},1);

                            for j=1:numels{1,1}
                                el = textscan(fid,'%f',1);
                                P{cur_shot,1}(j)= el{1,1};
                            end
                        elseif ~isempty(strfind(line, t_P_node))
                            datatype = textscan(fid,'%u',1);
                            numdims = textscan(fid,'%u',1);
                            numels = textscan(fid,'%u',1);

                            t_P{cur_shot,1}=zeros(numels{1,1},1);

                            for j=1:numels{1,1}
                                el = textscan(fid,'%f',1);
                                t_P{cur_shot,1}(j)= el{1,1};
                            end
                        end



                        line = fgets(fid);
                    end

                    fclose(fid);

                end
            end
        end
    end

    rawdata = struct('wl', wl, 'spec', I, 'int', int, 't_spec', t_I, 'P', P, 't_P', t_P, 'r', r);
end