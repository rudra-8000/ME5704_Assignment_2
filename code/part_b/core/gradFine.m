function [Tx, Ty] = gradFine(T, n, s)
% GRADFINE  Central-difference gradient (Table 5, O(h^2)) at the 9 given node positions
%   x, y = 1, 2, 3 on a refined field T = solveFine(n), grid spacing h_n = L/n, step s*h_n.
%   T(i,j) is at x = (i-1)*h_n, y = (j-1)*h_n (as returned by solveFine).
P = plateProblem();  hn = P.L/n;
if nargin < 3, s = 1; end
assert(s*hn < abs(P.xLine - round(P.xLine)), 'gradFine: stencil would reach the source line');
Tx = zeros(3);  Ty = zeros(3);
for j = 1:3, for i = 1:3
    I = round(i*P.h/hn) + 1;  J = round(j*P.h/hn) + 1;
    Tx(i,j) = (T(I+s,J) - T(I-s,J)) / (2*s*hn);
    Ty(i,j) = (T(I,J+s) - T(I,J-s)) / (2*s*hn);
end, end
end
