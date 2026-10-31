function P = plateProblem()
% PLATEPROBLEM  Problem data shared by parts (a), (b), (c).
%
%   Square plate [0,L]x[0,L], unit thickness, Laplace's equation inside.
%   All four edges held at Tb. A heat-source line (temperature Ts) lies on
%   x = xLine, y in [yLo, yHi] (closed segment, endpoints included).
%
%   The 9 interior grid nodes are (i*h, j*h), i,j = 1..3, with h = 1.
%   Convention: T(i,j) is the temperature at x = i*h, y = j*h  (i = x-index).

P.L     = 4;      % plate side length
P.h     = 1;      % grid spacing of the given data
P.Tb    = 200;    % edge temperature [K]
P.Ts    = 350;    % source-line temperature [K]
P.xLine = 1.5;    % source line x-position
P.yLo   = 2;      % source line lower end
P.yHi   = 3;      % source line upper end
P.nInt  = 3;      % interior nodes per side (given grid)
end
