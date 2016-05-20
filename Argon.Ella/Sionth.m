function Sionth = Sionth(Eion)
% n is not used => no need to have it as input
%function Sionth(n,Eion)
% thermalized Ar atoms: Tg=500 K -> E=3/2kTg=0.06 eV
% CHANGED to 0.03 eV AMK

%parameter(E=0.03)
%!parameter(E=0.00)
E=0.03;

Sionth=0.0;

%if(E.gt.Eion)then
if E > Eion
    bn=8.69e-18 * Eion^(-2.26);
    Sionth=bn * (E-Eion);        
end

%return
%end
