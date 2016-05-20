function Siona = Siona(E,n,Eion)

Siona=0.0;

%if(E.ge.Eion)then
% cfr vroeger (fit from exp)
if E >= Eion
    %if(n.eq.1)then
    if n == 1
    
        %if(E.le.75)Siona=10**(-29.175+6.554*log10(E))
        if E <= 75
            Siona=10^(-29.175 + 6.554*log10(E));
        end
        
        %if((E.gt.75).and.(E.le.100))Siona=10**(-23.7+3.636*log10(E))
        if (E > 75) & (E <= 100)
            Siona=10^(-23.7 + 3.636*log10(E));
        end
        
        %if((E.gt.100).and.(E.le.133.4))
        %.        Siona=10**(-20.125+1.8468*log10(E))
        if (E > 100) & (E <= 133.4)
            Siona=10^(-20.125 + 1.8468*log10(E));
        end
        
        %if((E.gt.133.4).and.(E.le.237))
        %.        Siona=10**(-18.518+1.0938*log10(E))
        if (E > 133.4) & (E <= 237)
            Siona=10^(-18.518 + 1.0938*log10(E));
        end
        
        %if((E.gt.237).and.(E.le.1000))
        %.        Siona=10**(-16.875+0.4018*log10(E))
        if (E > 237) & (E <= 1000)
            Siona=10^(-16.875 + 0.4018*log10(E));
        end
        
        %if(E.gt.1000)Siona=10**(-16.373+0.2346*log10(E))
        if E > 1000
            Siona=10^(-16.373 + 0.2346*log10(E));
        end
        
    else % n ~= 1
    % ion+atom ionization from level n (from Vlcek)
        xmea=2.725e-5;
        fac=7.3258e-17;
        Siona=fac/Eion^2 * (E/Eion-1)/(1+xmea*(E/Eion-1))^2;          
    end        
end

%return
%end
