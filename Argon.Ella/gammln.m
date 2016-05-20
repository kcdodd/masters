function gln = gammln(xx)
%function gammln(xx)
%gamma function added for error function 9/22/03 AMK	

%real gammln,xx
%INTEGER j
%DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)

% might need to have cof and std as outputs of gammln if they are used in
% another subroutine
%SAVE cof,stp

%DATA cof,stp/76.18009172947146D0,-86.50532032941677D0,
%     .	24.01409824083091D0,-1.231739572450155D0,.1208650973866179D-2,
%     .  -.5395239384953D-5,2.5066282746310005D0/
cof=[76.18009172947146e0 -86.50532032941677e0 24.01409824083091e0 -1.231739572450155e0 .1208650973866179e-2 -.5395239384953e-5];
stp=2.5066282746310005e0;

x=xx;
y=x;
tmp=x + 5.5e0;
tmp=(x + 0.5e0) * log(tmp) - tmp;
ser=1.000000000190015e0;

for j=1:6
    y=y + 1.e0;
    ser=ser + cof(j)/y;
end	

gln=tmp + log(stp * ser/x);

%return
%END