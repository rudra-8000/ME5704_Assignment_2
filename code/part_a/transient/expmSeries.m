function E = expmSeries(M)
% EXPMSERIES  Matrix exponential exp(M) by scaling and squaring with a Taylor series.
%   exp(M) = ( exp(M/2^s) )^(2^s); M/2^s is made small (norm < 0.5) so 20 Taylor terms
%   are accurate to machine precision.  Used as the exact reference for time marching:
%   for dT/dt = A*T - b the solution is  T(t) = Tss + exp(A*t) (T0 - Tss),  Tss = A\b.
n = size(M,1);  s = max(0, ceil(log2(max(1, max(sum(abs(M), 2)) / 0.5))));
X = M / 2^s;  E = eye(n);  term = eye(n);
for k = 1:20
    term = term * X / k;  E = E + term;
end
for k = 1:s, E = E * E; end
end
