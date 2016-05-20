function Sione_A = Sione(E_A,n,Eion)
% ionization of Ar0 (From Carman, J.Phys.D, 22, 55 ('89)):
% (= Bretagne et al., J.Phys.D, 14, 1225 ('81)):
% ionization of different Ar* levels (From Vlcek, J.Phys.D,22, 623 ('89))

N_E = size(E_A, 2);

Sione_A = zeros(1, N_E);

for E_i = 1:N_E
    
    E = E_A(E_i);

    Sione=0.0;

    %if(E.ge.Eion)then
    if E >= Eion
        %if(n.eq.1)then
        if n == 1
            a=(E-Eion)/2.0;
            EE=1.2 - 250.0/(E + 2*Eion);

            Sione=1e-16 * (23.9/E) * log((E+150.0/E)/Eion) * 4.6 * (atan((a-EE)/4.6) - atan(-EE/4.6));

        else
            %if((n.ge.2).and.(n.le.5))a=0.35
            if (n >= 2) && (n <= 5)
                a=0.35;
            end

            %if(n.eq.6)a=0.45
            if n == 6
                a=0.45;
            end

            %if((n.ge.7).and.(n.le.9))a=0.39
            if (n >= 7) && (n <= 9)
                a=0.39;
            end

            %if((n.ge.10).and.(n.le.11))a=0.32
            if (n >= 10) && (n <= 11)
                a=0.32;
            end

            %if(n.gt.11)a=0.67
            if n > 11
                a=0.67;
            end

            %if((n.ge.2).and.(n.le.11))b=4.0
            if (n >= 2) && (n <= 11)
                b=4.0;
            end

            %if(n.gt.11)b=1.0
            if n > 11
                b=1.0;
            end

            %Sione=6.783e-14/Eion**2*a*(Eion/E)**2*(E/Eion-1)*
            %.      log(1.25*b*E/Eion) 
            Sione=(6.783e-14/Eion^2) * a * (Eion/E)^2 * (E/Eion-1) * log(1.25 * b * E/Eion);        

        end        
    end
    
    Sione_A(E_i) = Sione;
    
end %E_i

%return
%end