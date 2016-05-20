function gammcf = gamcf(a,x)
%	gamma function added for error function 9/22/03 AMK
%	SUBROUTINE gcf(gammcf,a,x,gln)

%	INTEGER ITMAX
% no need to define the variable ITMAX=zeros since we assign a value to it
% later
%   REAL a,gammcf,gln,x,EPS,FPMIN
% we don't define a and x since they are given as inputs 
gammcf=zeros;
gln=zeros;
% no need to define the variable EPS=zeros and FPMIN=zeros since
% we assign values to them later
%	PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
ITMAX=100;
EPS=3e-7;
FPMIN=1e-30;

% 	USES gammln	
%	INTEGER i
% no need to define the variable i=zeros since we assign a value to it later
%	REAL an,b,c,d,del,h,gammln
% don't define gammln it's a function
% don't need to define the variable, we do so by assigning values to them
% later
gln=gammln(a);

b=x+1-a ;
c=1/FPMIN;
d=1/b;
h=d;

%do 912 i=1,ITMAX 
for i=1:ITMAX 
    an=-i * (i-a);
    b=b+2;
	d=an*d + b;
    %if(abs(d).lt.FPMIN)d=FPMIN
    if abs(d) < FPMIN
        d=FPMIN;
    end
    
    c=b + an/c;
    %if(abs(c).lt.FPMIN)c=FPMIN
    if abs(c) < FPMIN
        c=FPMIN;
    end
    
    d=1/d;
    del=d*c;
    h=h*del;
    
    %if(abs(del-1.).lt.EPS)goto 902
    if abs(del-1) < EPS
        break
    end
end

%pause 
%902 	gammcf=exp(-x+a*log(x)-gln)*h 
if abs(del-1) < EPS
    gammcf= exp(-x + a*log(x) - gln) * h;
else
    disp('error: gcf')
    pause 
end
%return
%END