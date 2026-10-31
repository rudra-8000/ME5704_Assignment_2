function [x, hist] = sor(A, b, omega, tol, maxIter)
% SOR  Successive over-relaxation (Liebmann's method for omega = 1 = Gauss-Seidel).
%   x_i^new = (1-omega)*x_i^old + omega * ( b_i - sum_{j~=i} A_ij x_j ) / A_ii
%   using already-updated values as soon as they exist.
%   hist(k) = max-norm residual after sweep k; stops when hist < tol.
n = numel(b);  x = zeros(n,1);  hist = zeros(maxIter,1);
for k = 1:maxIter
    for i = 1:n
        s = A(i,:)*x - A(i,i)*x(i);                      % sum over j ~= i with newest values
        x(i) = (1-omega)*x(i) + omega*(b(i) - s)/A(i,i);
    end
    hist(k) = max(abs(b - A*x));
    if hist(k) < tol, hist = hist(1:k); return, end
end
end
