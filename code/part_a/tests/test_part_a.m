% TEST_PART_A  Validation of the part (a) scripts against known answers.
%   1. Lecture-notes exercise (2x2 interior, SOR omega = 1.1)
%   2. Grid convergence on an exact harmonic function  u = e^x sin y   (expect order 2)
%   3. Patch test of the line treatments on  u = x^2 - y^2 + x*y  (quadratics: M1 and M3 exact)
%   4. All solvers agree with MATLAB's backslash (used ONLY here, as an independent oracle)
%   5. Discrete maximum principle for every treatment;  Lagrange weights reproduce a quadratic
nFail = 0;
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');

% ---- 1. lecture exercise: unit square, edges 0 (bottom) 100 (top) 75 (left) 50 (right), h = 1/3
disp('Test 1: lecture exercise');
A1 = [4 -1 -1 0; -1 4 0 -1; -1 0 4 -1; 0 -1 -1 4];          % order T11, T21, T12, T22
b1 = [0+75; 0+50; 75+100; 50+100];
x1 = sor(A1, b1, 1.1, 1e-10, 500);
fprintf('  T11 %.2f  T12 %.2f  T21 %.2f  T22 %.2f   (notes: 46.9, 71.9, 40.6, 65.6)\n', x1(1), x1(3), x1(2), x1(4));
nFail = nFail + ~check('matches the notes to 0.05 K', max(abs(x1 - [46.9; 40.6; 71.9; 65.6])) < 0.05);

% ---- 2. order of accuracy on u = exp(x) sin(y), unit square
disp('Test 2: grid convergence (exact harmonic solution)');
u = @(x,y) exp(x).*sin(y);  ns = [4 8 16 32];  err = zeros(size(ns));
for k = 1:numel(ns)
    n = ns(k);  [X, Y] = meshgrid((0:n)/n, (0:n)/n);  U = u(X', Y');   % U(i,j) at x = (i-1)/n
    T = zeros(n+1);  fixed = false(n+1);  fixed([1 end],:) = true;  fixed(:,[1 end]) = true;
    T(fixed) = U(fixed);
    T = sorGrid(T, fixed, 2/(1+sin(pi/n)), 1e-12, 20000);
    err(k) = max(max(abs(T - U)));
end
order = log2(err(1:end-1) ./ err(2:end));
fprintf('  h = 1/%s  max error = %s\n  observed order = %s\n', mat2str(ns), mat2str(err, 3), mat2str(order, 3));
nFail = nFail + ~check('observed order is 2 (+-0.1)', all(abs(order(2:end) - 2) < 0.1));
dlmwrite_csv(fullfile(resDir, 'test_convergence.csv'), [1./ns', err'], 'h,max_error');

% ---- 3. patch test of the source-line treatments
disp('Test 3: patch test, exact quadratic harmonic data on edges and line');
uq = @(x,y) x.^2 - y.^2 + x.*y;  exact = zeros(9,1);
for r = 1:9, i = mod(r-1,3)+1;  j = (r-i)/3 + 1;  exact(r) = uq(i, j); end
for m = {'sw', 'ghost3', 'ghost1'}
    [A, b] = buildSystem(m{1}, uq);
    e = max(abs(gaussElim(A, b) - exact));
    fprintf('  %-7s max error %.2e\n', m{1}, e);
    if ~strcmp(m{1}, 'ghost1'), nFail = nFail + ~check([m{1} ' reproduces the quadratic (error < 1e-9)'], e < 1e-9); end
end

% ---- 4. solver agreement
disp('Test 4: solvers vs backslash oracle');
[A, b] = buildSystem('ghost3');  xo = A\b;
nFail = nFail + ~check('Gauss elimination', max(abs(gaussElim(A, b) - xo)) < 1e-9);
nFail = nFail + ~check('Crout LU',          max(abs(luSolve(A, b) - xo)) < 1e-9);
nFail = nFail + ~check('Jacobi',            max(abs(jacobi(A, b, 1e-11, 5000) - xo)) < 1e-8);
nFail = nFail + ~check('Gauss-Seidel',      max(abs(sor(A, b, 1.0, 1e-11, 5000) - xo)) < 1e-8);
nFail = nFail + ~check('SOR 1.3',           max(abs(sor(A, b, 1.3, 1e-11, 5000) - xo)) < 1e-8);

% ---- 5. physics / interpolation sanity
disp('Test 5: maximum principle and Lagrange weights');
P = plateProblem();  ok = true;
for m = {'sw', 'ghost1', 'ghost3'}
    [A, b] = buildSystem(m{1});  x = gaussElim(A, b);  ok = ok && all(x >= P.Tb) && all(x <= P.Ts);
end
nFail = nFail + ~check('200 K <= T <= 350 K for SW and ghost treatments', ok);
w = lagrangeWeights([1 1.5 0], 2);
nFail = nFail + ~check('Lagrange weights reproduce f = x^2 and f = 1', abs(w*([1 1.5 0].^2)' - 4) < 1e-12 && abs(sum(w) - 1) < 1e-12);

fprintf('\n%d test group(s) failed\n', nFail);
if nFail > 0, error('part (a) tests failed'); end

function ok = check(name, cond)
% Print one PASS/FAIL line and return the result.
label = 'FAIL';  if cond, label = 'PASS'; end
fprintf('  [%s] %s\n', label, name);  ok = cond;
end

function dlmwrite_csv(file, M, header)
% Tiny CSV writer (header line + numeric matrix).
fid = fopen(file, 'w');  fprintf(fid, '%s\n', header);
for r = 1:size(M,1), fprintf(fid, '%.6g,', M(r,1:end-1));  fprintf(fid, '%.6g\n', M(r,end)); end
fclose(fid);
end
