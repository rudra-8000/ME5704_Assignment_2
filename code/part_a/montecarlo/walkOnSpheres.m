function [Test, se] = walkOnSpheres(x0, y0, nWalks, epsAbs, uEdge)
% WALKONSPHERES  Mesh-free Monte Carlo solution of Laplace's equation at (X0,Y0).
%   T(p) is harmonic, so it equals its average over any circle inside the domain.
%   Each walker jumps to a random point on the largest circle that touches the nearest
%   Dirichlet set (plate edge or source line) and repeats until it is within EPSABS of
%   one; it then scores that set's temperature.  T(x0,y0) = E[score].
%   No grid is used: the source line is treated exactly (distance to a segment).
%   Remaining errors: statistical (SE/sqrt(N)) and the absorption layer EPSABS.
%   Optional UEDGE = @(x,y): known edge values and no source line (used by the tests only).
P = plateProblem();
useLine = nargin < 5;
x = x0*ones(nWalks,1);  y = y0*ones(nWalks,1);  score = zeros(nWalks,1);  alive = true(nWalks,1);
while any(alive)
    xa = x(alive);  ya = y(alive);
    dEdge = min(min(xa, P.L-xa), min(ya, P.L-ya));                    % distance to plate edges
    dx    = abs(xa - P.xLine);  dy = max(max(P.yLo - ya, ya - P.yHi), 0);
    dLine = sqrt(dx.^2 + dy.^2);                                      % distance to the segment
    if ~useLine, dLine = Inf(size(dLine)); end
    d = min(dEdge, dLine);
    hit = d < epsAbs;                                                 % close enough: absorb
    if useLine
        sc = P.Tb*ones(size(d));  sc(dLine < dEdge) = P.Ts;
    else
        sc = uEdge(xa, ya);
    end
    idx = find(alive);
    score(idx(hit)) = sc(hit);
    th = 2*pi*rand(size(d));                                          % random jump on the circle
    x(idx) = xa + d.*cos(th);  y(idx) = ya + d.*sin(th);
    alive(idx(hit)) = false;
end
Test = mean(score);  se = std(score) / sqrt(nWalks);
end
