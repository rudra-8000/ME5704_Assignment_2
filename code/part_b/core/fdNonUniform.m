function w = fdNonUniform(a, c)
% FDNONUNIFORM  Three-point first derivative with unequal spacing a (left) and c (right):
%   f'(x) ~ w(1) f(x-a) + w(2) f(x) + w(3) f(x+c),   error -(a c / 6) f'''(x)  -> O(h^2).
%   Derived as in Lecture Notes Sec. 5.5.1: Taylor series at x-a and x+c, eliminate f''.
%   a = c gives the central difference C2 of Table 5.
w = [-c/(a*(a+c)), (c-a)/(a*c), a/(c*(a+c))];
end
