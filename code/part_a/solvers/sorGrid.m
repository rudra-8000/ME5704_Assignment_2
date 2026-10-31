function [T, nIter] = sorGrid(T, fixed, omega, tol, maxIter)
% SORGRID  Matrix-free SOR for Laplace's equation on a uniform grid.
%   T(i,j): temperatures, FIXED(i,j) = true marks Dirichlet nodes (never updated).
%   Each free node is relaxed towards the average of its 4 neighbours:
%       T_new = (1-omega)*T + omega*( T(i-1,j)+T(i+1,j)+T(i,j-1)+T(i,j+1) )/4
%   Stops when the largest update in a sweep is below TOL.
[nx, ny] = size(T);
for nIter = 1:maxIter
    maxChange = 0;
    for j = 2:ny-1
        for i = 2:nx-1
            if fixed(i,j), continue, end
            avg    = (T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1)) / 4;
            dT     = omega * (avg - T(i,j));
            T(i,j) = T(i,j) + dT;
            maxChange = max(maxChange, abs(dT));
        end
    end
    if maxChange < tol, return, end
end
end
