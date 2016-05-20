function Sdeexth = Sdeexth(m,n,Eexcm,Eexcn,gm,gn,gim,gin,g0m)

% thermalized Ar atoms: Tg=500 K -> E=3/2kTg=0.06 eV
% ion+atom deexcitation from level m to level n (from excitation)
% only between primed-primed, or unprimed-unprimed, no intercombination
% changed to 0.03; added missing Emn**(-2.26) 6/22/05 AMK
% modified to include 2-4,2-5,3-4,3-5 previously excluded by no intercomb AMK
% added abs to Emn AMK

%	parameter(E=0.00)
%   parameter(E=0.03)
E=0.03;
Sdeexth=0.0;
bmn=0.0;

%!if((E.gt.0.0).and.(n.lt.m))then
%!Emn=abs(Eexcm-Eexcn)
Emn=abs(Eexcm-Eexcn);

%if((E.gt.Emn).and.(n.lt.m))then
if (E > Emn) & (n < m)
%    Eext=E+Emn; NOT USED
    
    %if(((n.eq.2).and.((m.eq.4).or.(m.eq.5))).or.((n.eq.3)
    % .    .and.((m.eq.4).or.(m.eq.5))))bmn=4.8e-22*Emn**(-2.26)
    if ((n == 2) & ((m == 4) | (m == 5))) | ((n == 3) & ((m == 4) | (m == 5)))
        bmn=4.8e-22 * Emn^(-2.26);
    end
    
    %if(gin.eq.gim)then
    if gin == gim
        %if(((n.eq.2).and.(m.eq.3)).or.((n.eq.4).and.(m.eq.5)))then
        if ((n == 2) & (m == 3)) | ((n == 4) & (m == 5))
            bmn=1.79e-20 * Emn^(-2.26);
        %else if(((n.eq.2).and.((m.eq.4).or.(m.eq.5))).or.((n.eq.3)
        %.    .and.((m.eq.4).or.(m.eq.5))))then
        elseif ((n == 2) & ((m == 4) | (m == 5))) | ((n == 3) & ((m == 4) | (m == 5)))
            bmn=4.8e-22 * Emn^(-2.26);
        else
            bmn=8.69e-18 * Emn^(-2.26);
        end
    end

    %!Sext=bmn*(Eext-Emn)
    Sext=bmn * (E-Emn);
    
    %if(n.eq.1)gn=g0m
    if n == 1
        gn=g0m;
    end
    
    %!Sdeexth=gn/gm*Eext/E*Sext
    Sdeexth=gn/gm * E/(E-Emn) * Sext; 
end
%return
%end