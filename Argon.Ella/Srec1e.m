function Srec1e_A = Srec1e(E_A,n,Eion,gam,gn,gin)
% radiative recombination (from photoionization; inverse process)

% KCD vectorized E and return.

N_E = size(E_A, 2);

Srec1e_A = zeros(1, N_E);

for E_i = 1:N_E
    
    E = E_A(E_i);

    %parameter(E1h=13.884)
    E1h=13.884;
    Srec1e=0.0;
    hv=E+Eion;

    %if(E.gt.0.0)then
    if E > 0.0
        %if(n.eq.1)then
        if n == 1
            %if((hv.ge.Eion).and.(hv.le.(2*E1h)))Sp=3.5e-17
            if (hv >= Eion) && (hv <= (2*E1h))
                Sp=3.5e-17;
            end

            %if(hv.gt.(2*E1h))Sp=2.8e-16*(E1h/hv)**3
            if hv > (2*E1h)
                Sp=2.8e-16 * (E1h/hv)^3;
            end

        %else if((n.ge.2).and.(n.le.5))then
        elseif (n >= 2) && (n <= 5)
            %if((hv.ge.Eion).and.(hv.le.(0.59*E1h)))Sp=2e-18*gam
            if (hv >= Eion) && (hv <= (0.59*E1h))
                Sp=2e-18 * gam;
            end

            %if(hv.gt.(0.59*E1h))Sp=7.91e-18*gam*(Eion/E1h)**2.5*
            %.    (E1h/hv)**3
            if hv > (0.59*E1h)
                Sp=7.91e-18 * gam * (Eion/E1h)^2.5 * (E1h/hv)^3;
            end

        else % n > 5
            %if(hv.gt.Eion)Sp=gam*7.91e-18*(Eion/E1h)**2.5*(E1h/hv)**3
            if hv > Eion
                Sp=7.91e-18 * gam * (Eion/E1h)^2.5 * (E1h/hv)^3;
            end
        end

        Srec1e=(gn/(2*gin*5.1173e5)) * (hv^2/E) * Sp;        

    end
    
    Srec1e_A(E_i) = Srec1e;
    
end% E_i

%return
%end
