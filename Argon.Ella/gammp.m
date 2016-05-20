function gammp = gammp(a,x)
%gamma function added for error function 9/22/03 AMK

%REAL a,gammp,x
gammp=zeros;

%USES gcf,gser
%REAL gammcf,gamser,gln

%if(x.lt.0..or.a.le.0.)pause 
%if (x < 0) | (a <= 0) 
    %disp('gammp line 9')
    %x
    %a
%end

%if(x.lt.a+1.)then
if x < a+1
    %call gser(gamser,a,x,gln)
    gamser=gser(a,x);
    gammp=gamser;
else
    %call gcf(gammcf,a,x,gln)
    gammcf=gamcf(a,x);
    gammp=1-gammcf;
end 

%return
%end
