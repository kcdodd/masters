function Sexce_A = Sexce(E_A,n,m,Eexcn,Eexcm,aAf,aP,gn,gm)
% excitation of different Ar* levels to different Ar* levels (From Vlcek):

N_E = size(E_A, 2);

Sexce_A = zeros(1, N_E);

for E_i = 1:N_E
    
    E = E_A(E_i);

    Emn=Eexcm-Eexcn;
    Sexce=0.0;

    %if(E.ge.Emn)then
    if E >= Emn
        % CASES n-m = 2-3, 2-4, 2-5, 3-4, 3-5
        %if((n.ge.2).and.(n.le.3).and.(m.ge.3).and.(m.le.5))then
        if (n >= 2) & (n <= 3) & (m >= 3) & (m <= 5)
            %if((n.eq.2).and.(m.eq.3))Q=1.0
            if (n == 2) & (m == 3)
                Q=1.0;
            end
            if (m == 4) | (m == 5) % all other cases
                Q=0.1;
            end

            Sexce=gm/gn * (E-Emn)/E * 5.797e-15 * Q * (E-Emn)^(-0.54);

        % CASE n-m = 4-5
        %else if((n.eq.4).and.(m.eq.5))then
        elseif (n == 4) & (m == 5)

            Sexce=gm/gn * (E-Emn)/E * 8.111e-16 * (E-Emn)^(-1.04);

        %else if((aAf.eq.0.0).and.(aP.eq.0.0))then
        elseif (aAf == 0.0) & (aP == 0.0)

            Sexce=0.0;

        %ALL OTHER CASES
        else 
            %if(n.ge.1)b=1.0
            if n >= 1
                b=1.0;
            end
            %if((n.eq.1).and.((m.eq.3).or.(m.eq.5).or.(m.eq.15).or.
            %.    (m.eq.16)))b=4.0
            if (n == 1) & ((m == 3) | (m == 5) | (m == 15) | (m == 16))
                b=4.0;
            end
            %if((n.eq.1).and.(m.eq.17))b=2.0
            if (n == 1) & (m == 17)
                b=2.0;
            end
            %if((n.eq.1).and.((m.eq.20).or.(m.eq.21).or.(m.eq.26).or.
            %.    (m.eq.27).or.(m.eq.33)))b=1.0
            if (n == 1) & ((m == 20) | (m == 21) | (m == 26) | (m == 27) | (m == 33))
                b=1.0;
            end
            SexcA=6.783e-14/Emn^2 * aAf * (Emn/E)^2 * (E/Emn-1) * log(1.25*b*E/Emn);
            SexcF=3.519e-16 * aP * (Emn/E) * (1-(Emn/E));
            %if((n.eq.1).and.((m.eq.2).or.(m.eq.4).or.(m.eq.12).or.
            %.    (m.eq.13)))SexcF=3.519e-16*aP*(Emn/E)**3*(1-(Emn/E)**2)
            if (n == 1) & ((m == 2) | (m == 4) | (m == 12) | (m == 13))
                SexcF=3.519e-16 * aP * (Emn/E)^3 * (1-(Emn/E)^2);
            end

            Sexce=SexcA+SexcF; 

        end
    end
    
    Sexce_A(E_i) = Sexce;
end %E_i

%return
%end
