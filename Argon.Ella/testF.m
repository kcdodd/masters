function out = testF (in)
N=size(in,2);

out = zeros(1,N);
for in_i =  1:N;
    out(in_i) = sin(in(in_i));
end