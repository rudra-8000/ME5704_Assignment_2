function [A, b] = buildSystem(method, u)
% BUILDSYSTEM  Discrete Laplace system A*T = b on the given 3x3 interior grid (h = 1).
%
%   Rows are the finite-difference form of  T_xx + T_yy = 0  (units 1/length^2),
%   unknowns ordered r = (j-1)*3 + i  for the node (x,y) = (i,j).
%
%   The source line x = 1.5 (y in [2,3]) lies BETWEEN nodes (1,j) and (2,j), j = 2,3.
%   METHOD chooses how the line enters those four equations:
%     'sw'      Shortley-Weller: non-uniform 3-point T_xx from a Taylor expansion,
%               spacing 0.5 towards the line (value Ts there).
%     'ghost1'  Temporary ghost node on the far side of the line, filled by a linear
%               Lagrange interpolation through (node, line).
%     'ghost3'  Same ghost node, cubic Lagrange through (node, line, opposite neighbour, and the
%               point two spacings beyond the node towards the line).
%     'snapL'   Naive: clamp nodes (1,2),(1,3) to Ts.
%     'snapR'   Naive: clamp nodes (2,2),(2,3) to Ts.
%     'snapAll' Naive: clamp all four nodes to Ts.
%   Optional U = @(x,y) exact harmonic function: edge and line values are then taken from it
%   (used only by the tests); by default edges = Tb and line = Ts.
P = plateProblem();  n = P.nInt;
edge = @(i,j) P.Tb;  line = @(j) P.Ts;
if nargin > 1, edge = @(i,j) u(i*P.h, j*P.h);  line = @(j) u(P.xLine, j*P.h); end
A = zeros(n*n);  b = zeros(n*n, 1);
clamp = {};                                   % list of clamped (i,j) for the snap methods
if strcmp(method, 'snapL'),   clamp = {[1 2], [1 3]}; end
if strcmp(method, 'snapR'),   clamp = {[2 2], [2 3]}; end
if strcmp(method, 'snapAll'), clamp = {[1 2], [1 3], [2 2], [2 3]}; end

for j = 1:n
    for i = 1:n
        r = (j-1)*n + i;
        if any(cellfun(@(c) isequal(c, [i j]), clamp))     % clamped row: T = Ts
            A(r,r) = 1;  b(r) = line(j);  continue
        end
        % --- y-direction: ordinary central difference (the line never crosses it)
        [A,b] = addT(A,b,r,i,j-1, 1, edge);  [A,b] = addT(A,b,r,i,j,-2, edge);  [A,b] = addT(A,b,r,i,j+1, 1, edge);
        % --- x-direction: is the line the nearest neighbour on the east (+1) or west (-1)?
        side = 0;
        if any(j == [2 3]) && i == 1, side = +1; end
        if any(j == [2 3]) && i == 2, side = -1; end
        if side == 0 || ~any(strcmp(method, {'sw','ghost1','ghost3'}))
            [A,b] = addT(A,b,r,i-1,j, 1, edge);  [A,b] = addT(A,b,r,i,j,-2, edge);  [A,b] = addT(A,b,r,i+1,j, 1, edge);
        elseif strcmp(method, 'sw')
            dm = 1;  dp = 1;                                 % distances to W and E neighbours
            if side > 0, dp = 0.5; else, dm = 0.5; end
            wm = 2/(dm*(dm+dp));  wp = 2/(dp*(dm+dp));       % T_xx ~ wm*Tm - (wm+wp)*T + wp*Tp
            if side > 0
                [A,b] = addT(A,b,r,i-1,j,wm, edge);  b(r) = b(r) - wp*line(j);
            else
                [A,b] = addT(A,b,r,i+1,j,wp, edge);  b(r) = b(r) - wm*line(j);
            end
            [A,b] = addT(A,b,r,i,j,-(wm+wp), edge);
        else                                                  % ghost1 / ghost3
            % ghost value at x = i + side, interpolated from the points xs:
            %   node i, line, opposite neighbour i-side, (cubic only) the point i+2*side
            if strcmp(method, 'ghost1'), xs = [i, P.xLine]; else, xs = [i, P.xLine, i-side, i+2*side]; end
            w = lagrangeWeights(xs, i + side);  w(end+1:4) = 0;   % ghost = w1*T_P + w2*line + w3*T_opp + w4*T_far
            [A,b] = addT(A,b,r,i-side,j,   1 + w(3), edge);       % opposite neighbour + its share of the ghost
            [A,b] = addT(A,b,r,i,j,       -2 + w(1), edge);
            [A,b] = addT(A,b,r,i+2*side,j,     w(4), edge);       % far point (zero weight for ghost1)
            b(r) = b(r) - w(2)*line(j);
        end
    end
end
end

function [A, b] = addT(A, b, r, i, j, c, edge)
% Add c*T(i,j) to row r; a node outside the interior is an edge node (T = Tb) -> moves to b.
P = plateProblem();  n = P.nInt;
if i >= 1 && i <= n && j >= 1 && j <= n
    A(r, (j-1)*n + i) = A(r, (j-1)*n + i) + c;
else
    b(r) = b(r) - c*edge(i,j);
end
end
