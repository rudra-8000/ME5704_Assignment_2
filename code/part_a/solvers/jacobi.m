function [x, hist] = jacobi(A, b, tol, maxIter)
% JACOBI  Jacobi iteration: every unknown is updated from the previous sweep only.
%   x_i^(k+1) = ( b_i - sum_{j~=i} A_ij x_j^(k) ) / A_ii
%   hist(k) = max-norm residual ||b - A x||_inf after sweep k; stops when hist < tol.
n = numel(b);  x = zeros(n,1);  hist = zeros(maxIter,1);  d = diag(A);
for k = 1:maxIter
    x = (b - (A*x - d.*x)) ./ d;
    hist(k) = max(abs(b - A*x));
    if hist(k) < tol, hist = hist(1:k); return, end
end
end
