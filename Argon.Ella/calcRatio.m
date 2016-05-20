function [n0_ratio, nion_ratio] = calcRatio(Te)

ArFracs = [
            1.5 1 0 0 0 0;
            1.7 0.999 0.001 0 0 0;
            2.0 0.973 0.027 0 0 0;
            2.5 0.801 0.199 0 0 0;
            3.0 0.629 0.368 0.003 0 0;
            4.0 0.385 0.551 0.064 0 0;
            4.5 0.302 0.568 0.130 0 0;
            5.0 0.213 0.543 0.244 0 0;
            6.1 0.110 0.431 0.458 0.002 0;
            7.5 0.043 0.272 0.671 0.013 0;
            8.0 0.033 0.232 0.715 0.021 0;
            8.9 0.019 0.168 0.768 0.046 0;
            9.8 0.012 0.126 0.782 0.080 0.001;
            10.5 0.008 0.100 0.775 0.116 0.002;
            11.4 0.005 0.074 0.743 0.175 0.004];
        
for t = 1:15
   if (Te >= ArFracs(t,1) && Te <= ArFracs(t+1,1))
       diff = (ArFracs(t+1,1)-ArFracs(t,1));
       frac = (Te - ArFracs(t,1))/diff;
       
       ch_sum = 0;
       
       for z = 1:5
           ch_sum = ch_sum + (z-1)*((1-frac)*ArFracs(t,z+1) + frac*ArFracs(t+1,z+1));
       end
       
       n0_ratio = ((1-frac)*ArFracs(t,2) + frac*ArFracs(t+1,2))/ch_sum;
       nion_ratio = (1 - ((1-frac)*ArFracs(t,2) + frac*ArFracs(t+1,2)))/ch_sum;
       
       break;
   end
end