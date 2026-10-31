function w = lagrangeWeights(xs, xg)
% LAGRANGEWEIGHTS  Weights w so that p(xg) = sum_k w(k)*f(xs(k)), where p is
%   the Lagrange polynomial through the points (xs(k), f(xs(k))).
%
%       L_k(x) = prod_{m ~= k} (x - xs(m)) / (xs(k) - xs(m))
%
%   Used to fill a temporary "ghost" node by interpolation (degree = numel(xs)-1).
w = ones(size(xs));
for k = 1:numel(xs)
    for m = 1:numel(xs)
        if m ~= k
            w(k) = w(k) * (xg - xs(m)) / (xs(k) - xs(m));
        end
    end
end
end
