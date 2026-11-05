function D = gradGiven(T9, u)
% GRADGIVEN  Temperature gradient at the 9 given nodes (h = 1) by every scheme compared in part (b).
%   T9(i,j) = temperature at x = i*h, y = j*h (i,j = 1..3), e.g. column T_cubic_h1 of part (a).
%   Edge values are Tb and the line value is Ts, unless a function U(x,y) is given (used by the
%   tests): then edges and line take the values of U.
%   D.(scheme).Tx, D.(scheme).Ty are 3x3 arrays [K/m]; NaN where the scheme cannot be used.
%     C2   central O(h^2), Table 5. At the four nodes beside the line it crosses the line (= N0).
%     C4   central O(h^4), Table 5, where x +- 2h stays on one side of the line and on the plate.
%     R    Richardson 4/3 D(h) - 1/3 D(2h) with C2 (Sec. 5.5.2); algebraically equal to C4.
%     OS1  one-sided O(h) on the side away from the line (B1 left of it, F1 right of it).
%     OS2  one-sided O(h^2) away from the line (B2 left / F2 right); needs two nodes on that side.
%     NU   three points including the line value at distance |x - xLine| (fdNonUniform, O(h^2)).
P = plateProblem();  h = P.h;  n = P.L/h;               % n = 4 divisions per side
if nargin < 2, u = []; end
G = P.Tb * ones(n+1);                                   % G(i+1,j+1) = T at x = i*h, y = j*h
lineVal = @(y) P.Ts;
if ~isempty(u)
    [X, Y] = ndgrid((0:n)*h);  G = u(X, Y);  lineVal = @(y) u(P.xLine, y);
end
G(2:n, 2:n) = T9;

names = {'C2', 'C4', 'R', 'OS1', 'OS2', 'NU'};
for s = names, D.(s{1}).Tx = nan(3);  D.(s{1}).Ty = nan(3); end
for j = 1:3, for i = 1:3
    x = i*h;  y = j*h;  I = i+1;  J = j+1;
    onLineRow = y >= P.yLo && y <= P.yHi;
    % ---- y-derivatives: the line is vertical and between node columns, so no y-stencil crosses it
    D.C2.Ty(i,j) = applyFd(G(I,:), J, 'C2', h, 1);
    D.C4.Ty(i,j) = applyFd(G(I,:), J, 'C4', h, 1);
    D.R.Ty(i,j)  = richardson(G(I,:), J, h);
    % ---- x-derivatives
    D.C2.Tx(i,j) = applyFd(G(:,J), I, 'C2', h, 1);
    if ~(onLineRow && abs(x - P.xLine) < 2*h)          % +-2h stencil stays on one side of the line
        D.C4.Tx(i,j) = applyFd(G(:,J), I, 'C4', h, 1);
        D.R.Tx(i,j)  = richardson(G(:,J), I, h);
    end
    if onLineRow && abs(x - P.xLine) < h               % the +-h stencil crosses the line
        d = P.xLine - x;                               % signed distance from node to line
        if d > 0                                       % line to the right: use the left side
            D.OS1.Tx(i,j) = applyFd(G(:,J), I, 'B1', h, 1);
            D.OS2.Tx(i,j) = applyFd(G(:,J), I, 'B2', h, 1);    % NaN at x = h (would need x = -h)
            D.NU.Tx(i,j)  = fdNonUniform(h, d) * [G(I-1,J); G(I,J); lineVal(y)];
        else                                           % line to the left: use the right side
            D.OS1.Tx(i,j) = applyFd(G(:,J), I, 'F1', h, 1);
            D.OS2.Tx(i,j) = applyFd(G(:,J), I, 'F2', h, 1);
            D.NU.Tx(i,j)  = fdNonUniform(-d, h) * [lineVal(y); G(I,J); G(I+1,J)];
        end
    end
end, end
end

function d = applyFd(v, k, name, h, stride)
% Apply scheme NAME at index K of the vector V with step STRIDE*h; NaN if it leaves the plate.
[off, w] = fdWeights(name, stride*h);
idx = k + stride*off;
if any(idx < 1 | idx > numel(v)), d = NaN; return, end
d = w * reshape(v(idx), [], 1);
end

function d = richardson(v, k, h)
% Lecture Notes Sec. 5.5.2 with h1 = 2h, h2 = h:  D = 4/3 D(h) - 1/3 D(2h)  for an O(h^2) scheme.
d = 4/3 * applyFd(v, k, 'C2', h, 1) - 1/3 * applyFd(v, k, 'C2', h, 2);
end
