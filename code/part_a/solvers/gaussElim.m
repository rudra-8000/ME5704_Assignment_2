function x = gaussElim(A, b)
% GAUSSELIM  Solve A*x = b by Gauss elimination with partial pivoting.
%   Forward elimination to upper-triangular form, then back substitution.
n = numel(b);
for k = 1:n-1
    [~, p] = max(abs(A(k:n, k)));  p = p + k - 1;       % pivot row = largest |entry| in column k
    A([k p], :) = A([p k], :);  b([k p]) = b([p k]);
    for r = k+1:n
        f = A(r,k) / A(k,k);                            % multiplier
        A(r,k:n) = A(r,k:n) - f*A(k,k:n);
        b(r)     = b(r)     - f*b(k);
    end
end
x = zeros(n,1);
for k = n:-1:1                                          % back substitution
    x(k) = (b(k) - A(k,k+1:n)*x(k+1:n)) / A(k,k);
end
end
