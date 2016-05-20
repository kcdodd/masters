function Sreci = Sreci(E,n,Eion,gn,gin)
% ion+atom recomb to level n (from ionization:Siona)

% ! will be not used in MC
% !!!Srec=sigma/ne

%parameter(Ee=10.0,fac32=3.313e-22)
Ee=10.0;
fac32=3.313e-22;
Sreci=0.0;

%if(E.gt.0.0)then
if E > 0.0
    xkT=0.667 * Ee;
    fac32b=fac32 / xkT^1.5;
    E2=E + Eion;
    
    Sreci=gn/(2*gin) * fac32b * E2/E * Siona(E2,n,Eion);        
end

%return
%end
