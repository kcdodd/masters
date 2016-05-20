function Sdeexe_A = Sdeexe(E_A,m,n,Eexcm,Eexcn,aAf,aP,gm,gn,g0m)
% superelastic (deexcitation) collis.from Ar*(m) to Ar*(n): From excit 
% (detailed balancing):

N_E = size(E_A, 2);

Sdeexe_A = zeros(1, N_E);

for E_i = 1:N_E
    
    E = E_A(E_i);

    Sdeexe=0.0;

    %if(E.gt.0.0)then
    if E > 0.0
        Emn=Eexcm-Eexcn;
        %if(n.eq.1)gn=g0m
        if n == 1
            gn=g0m;
        end

        Sdeexe=(gn/gm) * ((E+Emn)/E) * Sexce(E+Emn,n,m,Eexcn,Eexcm,aAf,aP,gn,gm);

    end
    
    Sdeexe_A(E_i) = Sdeexe;
end %E_i
%return
%end
