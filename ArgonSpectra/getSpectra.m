function [wavelengths, avg] = getSpectra(tag, n1, n2)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 13-Oct-2011 15:44:17

N = n2 - n1 + 1;
% Import the file
for i = 1:N
    filename = [tag num2str(i + n1 - 1, '%05u') '.txt'];
    newData1 = importdata(filename);

    spectra = newData1.('data');
    
    if (i == 1)
       wavelengths = spectra(:,1);
       avg = spectra(:,2);
    else
        avg = avg + spectra(:,2);
    end

end

avg = avg/N;