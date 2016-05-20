function gamser = gser(a,x)
%function gser(gamser,a,x,gln)
%gamma function added for error function 9/22/03 AMK

%	INTEGER ITMAX
%	REAL a,gamser,gln,x,EPS
gamser = zeros;
gln=zeros;

%	PARAMETER (ITMAX=100,EPS=3.e-7)
ITMAX=100;
EPS=3.e-7;

%USES gammln
%	INTEGER n
%	REAL ap,del,sum,gammln

gln=gammln(a);

%if(x.le.0.)then
if x <= 0
    %if(x.lt.0.)pause 
    %if x < 0
        %disp('gser: x<0')
    %end
    gamser=0;
    return %exit the function
end

% if x>0 then we continue
ap=a;
sum=1/a;
del=sum;

%do 911 n=1,ITMAX
for n=1:ITMAX
    ap=ap + 1;
    del=del * x/ap;
    sum=sum + del;
    %if(abs(del).lt.abs(sum)*EPS)goto 901
    if abs(del) < abs(sum)*EPS
        break
    end
end
%pause
%901 	gamser=sum*exp(-x+a*log(x)-gln)
if abs(del) < abs(sum)*EPS
    gamser=sum * exp(-x + a * log(x) - gln);
else
    disp('gser:error')
    pause
end

%return
%end