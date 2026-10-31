function [T, nIter] = solveFine(n, omega, yLo, yHi)
% SOLVEFINE  Plate solution on a uniform grid with N divisions per side (h = L/N).
%   For N = 8, 16, 32, ... the source line x = 1.5, y in [2,3] falls exactly on grid
%   nodes, so it is imposed as ordinary Dirichlet nodes (no special stencil needed).
%   OMEGA = SOR factor (default: the classic optimum for a square, 2/(1+sin(pi/N));
%   omega = 1 gives plain Gauss-Seidel).  YLO, YHI optionally move the line ends
%   (used only by the robustness study; they must be multiples of h).
%   Returns T(i,j) with x = (i-1)*h, y = (j-1)*h.
P = plateProblem();  h = P.L / n;
if nargin < 3, yLo = P.yLo;  yHi = P.yHi; end
T = P.Tb * ones(n+1);                               % edges and initial guess = Tb
fixed = false(n+1);  fixed([1 end], :) = true;  fixed(:, [1 end]) = true;
iL = round(P.xLine/h) + 1;                          % column of the line
jL = round(yLo/h) + 1 : round(yHi/h) + 1;           % rows covered by the line
T(iL, jL) = P.Ts;  fixed(iL, jL) = true;
if nargin < 2 || isempty(omega), omega = 2 / (1 + sin(pi/n)); end
[T, nIter] = sorGrid(T, fixed, omega, 1e-10, 50000);
end
