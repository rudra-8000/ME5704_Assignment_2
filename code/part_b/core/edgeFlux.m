function [Qk, perEdge] = edgeFlux(T, h, order, lineSafe)
% EDGEFLUX  Heat leaving through the four edges per unit conductivity and unit thickness:
%   Q/k = sum over the edges of  integral of (-dT/dn) ds     [K m for T in K and lengths in m].
%   T(i,j) on a uniform grid with spacing h, x = (i-1)h, y = (j-1)h, edge nodes included.
%   dT/dn (outward normal) at every edge node, corners included, from the edge node and the
%   nodes inward along the normal (distances 0, h, 2h):
%     ORDER = 2 (default): 3-point one-sided O(h^2), Tables 3/4:  dT/dn = (3T_0 - 4T_1 + T_2)/(2h)
%     ORDER = 1          : 2-point one-sided O(h),   Tables 3/4:  dT/dn = (T_0 - T_1)/h
%   LINESAFE (default true): where the 3-point stencil would cross the source line, use the
%   2-point formula, which stays on the edge side of the line.
%   Along each edge the composite trapezoidal rule (Lecture Notes Sec. 5.2) integrates -dT/dn.
%   perEdge = [east west north south].
if nargin < 3, order = 2; end
if nargin < 4, lineSafe = true; end
P = plateProblem();  n = size(T,1) - 1;  y = (0:n)*h;
% profiles along the inward normal: rows = edge nodes, columns = distance 0, h, 2h
prof = {[T(end,:)', T(end-1,:)', T(end-2,:)'], ...      % east  (x = L)
        [T(1,:)',   T(2,:)',     T(3,:)'], ...          % west  (x = 0)
        [T(:,end),  T(:,end-1),  T(:,end-2)], ...       % north (y = L)
        [T(:,1),    T(:,2),      T(:,3)]};              % south (y = 0)
xEdge = [P.L, 0];                                       % x of the east and west edges
perEdge = zeros(1, 4);
for e = 1:4
    V = prof{e};
    dTdn1 = (V(:,1) - V(:,2)) / h;
    dTdn2 = (3*V(:,1) - 4*V(:,2) + V(:,3)) / (2*h);
    dTdn  = dTdn2;  if order == 1, dTdn = dTdn1; end
    if order == 2 && lineSafe && e <= 2                 % only x-stencils can cross the vertical line
        dist  = abs(P.xLine - xEdge(e));                % distance from this edge to the line
        cross = dist < 2*h & y' >= P.yLo & y' <= P.yHi;
        assert(~any(cross & dist <= h), 'edgeFlux: grid too coarse, 2-point stencil reaches the line');
        dTdn(cross) = dTdn1(cross);
    end
    perEdge(e) = trapezoid(-dTdn, h);
end
Qk = sum(perEdge);
end

function I = trapezoid(f, h)
% Composite trapezoidal rule, Lecture Notes Sec. 5.2:  h/2 [f_0 + 2 sum f_i + f_N].
I = h * (sum(f) - (f(1) + f(end))/2);
end
