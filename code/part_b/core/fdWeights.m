function [off, w] = fdWeights(name, h)
% FDWEIGHTS  First-derivative finite-difference schemes of Lecture Notes Tables 3-5 (Sec. 5.5.1).
%   f'(x) ~ sum_k w(k) * f(x + off(k)*h).
%   Errors: F1, B1 O(h);  C2, F2, B2 O(h^2);  C4 O(h^4).
switch name
    case 'F1', off = [0 1];        w = [-1 1];
    case 'B1', off = [-1 0];       w = [-1 1];
    case 'F2', off = [0 1 2];      w = [-3 4 -1] / 2;       % Table 3: [-f(i+2) + 4f(i+1) - 3f(i)] / 2h
    case 'B2', off = [-2 -1 0];    w = [1 -4 3] / 2;        % Table 4: [3f(i) - 4f(i-1) + f(i-2)] / 2h
    case 'C2', off = [-1 1];       w = [-1 1] / 2;          % Table 5: [f(i+1) - f(i-1)] / 2h
    case 'C4', off = [-2 -1 1 2];  w = [1 -8 8 -1] / 12;    % Table 5: [-f(i+2) + 8f(i+1) - 8f(i-1) + f(i-2)] / 12h
    otherwise, error('fdWeights: unknown scheme %s', name);
end
w = w / h;
end
