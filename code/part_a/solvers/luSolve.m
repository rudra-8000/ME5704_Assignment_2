function x = luSolve(A, b)
% LUSOLVE  Solve A*x = b by Crout LU decomposition (A = L*U, U has unit diagonal).
%   Step 1: factor A.  Step 2: L*y = b (forward).  Step 3: U*x = y (backward).
%   No pivoting: fine for the diagonally dominant Laplace matrices used here.
n = numel(b);  L = zeros(n);  U = eye(n);
for j = 1:n
    for i = j:n
        L(i,j) = A(i,j) - L(i,1:j-1)*U(1:j-1,j);
    end
    for i = j+1:n
        U(j,i) = (A(j,i) - L(j,1:j-1)*U(1:j-1,i)) / L(j,j);
    end
end
y = zeros(n,1);
for i = 1:n,    y(i) = (b(i) - L(i,1:i-1)*y(1:i-1)) / L(i,i); end
x = zeros(n,1);
for i = n:-1:1, x(i) = y(i) - U(i,i+1:n)*x(i+1:n); end
end
