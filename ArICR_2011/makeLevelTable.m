function [level_table] = makeLevelTable (n0, rates)

Te_l = [0.5 1 2 5 10 20];%eV
ne_l = 1e6*[1e10 2e10 5e10 1e11 2e11 5e11 1e12 2e12];%m^-3

N_Te = size(Te_l,2);
N_ne = size(ne_l,2);

% level 1=6, 2=7, 3=8, 4=9
npop = cell(4,1);

npop{1,1} = zeros(N_Te, N_ne);
npop{2,1} = zeros(N_Te, N_ne);
npop{3,1} = zeros(N_Te, N_ne);
npop{4,1} = zeros(N_Te, N_ne);

for Te_i = 1:N_Te
    for ne_i = 1:N_ne
        solution = arI_solution(Te_l(Te_i), ne_l(ne_i)/1e6, n0/1e6, 1, 2, rates, 0);
        
        npop{1,1}(Te_i,ne_i) = 1e6*solution.npop(6)/n0;
        npop{2,1}(Te_i,ne_i) = 1e6*solution.npop(7)/n0;
        npop{3,1}(Te_i,ne_i) = 1e6*solution.npop(8)/n0;
        npop{4,1}(Te_i,ne_i) = 1e6*solution.npop(9)/n0;
    end
    
end


level_table = struct('npop', npop, 'Te', Te_l, 'ne', ne_l);