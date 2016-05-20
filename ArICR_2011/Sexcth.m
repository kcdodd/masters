function Sexcth = Sexcth(n,m,Eexcn,Eexcm,gin,gim)
% thermalized Ar atoms: Tg=450 K -> E=3/2kTg=0.06 eV CHANGED TO 0.03eV
% only between primed-primed, or unprimed-unprimed, no intercombination
% modified to include 2-4,2-5,3-4,3-5 previously excluded by no intercomb AMK
% added abs to Emn AMK

%   parameter(E=0.03)
%!	parameter(E=0.00)
E=0.03;
Sexcth=0.0;
bmn=0.0;
Emn=abs(Eexcm-Eexcn);

%if((E.gt.Emn).and.(n.lt.m))then
if (E > Emn) & (n < m)
    %if(((n.eq.2).and.((m.eq.4).or.(m.eq.5))).or.((n.eq.3).and.
    % .  ((m.eq.4).or.(m.eq.5))))bmn=4.8e-22*Emn**(-2.26)
    if ((n == 2) & ((m == 4) | (m == 5))) | ((n == 3) & ((m == 4) | (m == 5)))
        bmn=4.8e-22 * Emn^(-2.26);
    end
    
    %if(gin.eq.gim)then
    if gin == gim
        %if(((n.eq.2).and.(m.eq.3)).or.((n.eq.4).and.(m.eq.5)))then
        if ((n == 2) & (m == 3)) | ((n == 4) & (m == 5))
            bmn=1.79e-20 * Emn^(-2.26);
        %else if(((n.eq.2).and.((m.eq.4).or.(m.eq.5))).or.((n.eq.3).and.
        %.  ((m.eq.4).or.(m.eq.5))))then
        elseif ((n == 2) & ((m == 4) | (m == 5))) | ((n == 3) & ((m == 4) | (m == 5)))
            bmn=4.8e-22 * Emn^(-2.26);
        else
            bmn=8.69e-18 * Emn^(-2.26);        
        end        
    end

    Sexcth=bmn * (E-Emn);        

end

%return
%end
