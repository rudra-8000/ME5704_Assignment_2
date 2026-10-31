% TEST_PART_A  Validation of the part (a) scripts against known answers.
%   1. Lecture-notes exercise (2x2 interior, SOR omega = 1.1)                  -> sor
%   2. Grid convergence on an exact harmonic function u = e^x sin y (order 2)  -> sorGrid
%   3. Patch tests of the line treatments: quadratic (M1, M3 exact) and cubic (only M3 exact) -> buildSystem
%   4. Slit-tip solution u = Re sqrt(z - z0): order 1 from the sqrt(r) singularity -> sorGrid on a slit
%   5. All solvers agree with MATLAB's backslash (used ONLY here, as an independent oracle)
%   6. Time-marching orders against the exact matrix-exponential solution      -> marchHeat
%   7. Monte Carlo against exact answers                                       -> walkGrid, walkOnSpheres
%   8. Discrete maximum principle;  Lagrange weights reproduce a quadratic
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

% ---- 3. patch tests of the source-line treatments (edge and line values taken from u)
disp('Test 3: patch tests, exact harmonic data on edges and line');
patch = {@(x,y) x.^2 - y.^2 + x.*y, 'quadratic x^2-y^2+xy'; @(x,y) x.^3 - 3*x.*y.^2, 'cubic x^3-3xy^2'};
patchErr = zeros(2, 3);
for q = 1:2
    exact = zeros(9,1);
    for r = 1:9, i = mod(r-1,3)+1;  j = (r-i)/3 + 1;  exact(r) = patch{q,1}(i, j); end
    for m = 1:3
        meth = {'sw', 'ghost1', 'ghost3'};
        [A, b] = buildSystem(meth{m}, patch{q,1});
        patchErr(q,m) = max(abs(gaussElim(A, b) - exact));
    end
    fprintf('  %-22s max error  M1 %.2e   M2 %.2e   M3 %.2e\n', patch{q,2}, patchErr(q,:));
end
nFail = nFail + ~check('quadratic: M1 and M3 exact (< 1e-9)', patchErr(1,1) < 1e-9 && patchErr(1,3) < 1e-9);
nFail = nFail + ~check('cubic: M3 exact (< 1e-9) while M1 and M2 are not', patchErr(2,3) < 1e-9 && patchErr(2,1) > 1e-3 && patchErr(2,2) > 1e-3);

% ---- 4. slit tip: u = Re sqrt(z - z0) is harmonic, zero on the slit y = 2, x <= 2, and behaves like sqrt(r) at the tip
disp('Test 4: slit-tip solution (sqrt(r) singularity)');
ns = [16 32 64 128];  errS = zeros(size(ns));
for k = 1:numel(ns)
    n = ns(k);  x = 4*(0:n)/n;  [X, Y] = meshgrid(x, x);  U = real(sqrt(complex(X' - 2, Y' - 2)));
    fixed = false(n+1);  fixed([1 end],:) = true;  fixed(:,[1 end]) = true;
    fixed(x <= 2, n/2 + 1) = true;                                % the slit nodes (y = 2, x <= 2)
    T = zeros(n+1);  T(fixed) = U(fixed);
    T = sorGrid(T, fixed, 2/(1+sin(pi/n)), 1e-11, 50000);
    pts = [3 2; 2 3; 3 3; 2 1; 1 3] * n/4 + 1;                    % fixed points at distance >= 1 from the tip
    errS(k) = max(abs(T(sub2ind(size(T), pts(:,1), pts(:,2))) - U(sub2ind(size(T), pts(:,1), pts(:,2)))));
end
orderS = log2(errS(1:end-1) ./ errS(2:end));
dlmwrite_csv(fullfile(resDir, 'test_slit.csv'), [4./ns', errS'], 'h,max_error');
fprintf('  max error at 5 fixed points = %s\n  observed order = %s\n', mat2str(errS, 3), mat2str(orderS, 3));
nFail = nFail + ~check('order is about 1 (0.7 to 1.3), not 2', all(orderS(2:end) > 0.7 & orderS(2:end) < 1.3));

% ---- 5. solver agreement
disp('Test 5: solvers vs backslash oracle');
[A, b] = buildSystem('ghost3');  xo = A\b;
nFail = nFail + ~check('Gauss elimination', max(abs(gaussElim(A, b) - xo)) < 1e-9);
nFail = nFail + ~check('Crout LU',          max(abs(luSolve(A, b) - xo)) < 1e-9);
nFail = nFail + ~check('Jacobi',            max(abs(jacobi(A, b, 1e-11, 5000) - xo)) < 1e-8);
nFail = nFail + ~check('Gauss-Seidel',      max(abs(sor(A, b, 1.0, 1e-11, 5000) - xo)) < 1e-8);
nFail = nFail + ~check('SOR 1.3',           max(abs(sor(A, b, 1.3, 1e-11, 5000) - xo)) < 1e-8);

% ---- 6. time marching: error at s = 1 for steps ds and ds/2 against Tss + exp(A s)(T0 - Tss)
disp('Test 6: time-marching orders vs matrix exponential');
Tss = gaussElim(A, b);  T0 = 200*ones(9,1);  sEnd = 1;
Texact = Tss + expmSeries(A*sEnd) * (T0 - Tss);
schemes = {'euler', 'ieuler', 'cn', 'rk4'};  expected = [1 1 2 4];  ratio = zeros(1, 4);
for q = 1:4
    e = zeros(1, 2);
    for k = 1:2
        ns = 10*2^k;  Th = marchHeat(A, b, T0, sEnd/ns, ns, schemes{q});   % ds = 0.05, 0.025
        e(k) = max(abs(Th(:,end) - Texact));
    end
    ratio(q) = e(1)/e(2);
    fprintf('  %-6s error %.2e -> %.2e, ratio %.1f (expected about %d)\n', schemes{q}, e, ratio(q), 2^expected(q));
    nFail = nFail + ~check([schemes{q} ' converges with order ' num2str(expected(q))], abs(log2(ratio(q)) - expected(q)) < 0.4);
end

% ---- 7. Monte Carlo against exact answers (quadratic harmonic function, so M3 and the continuum agree)
disp('Test 7: Monte Carlo vs exact harmonic function');
rng(2);  uq = patch{1,1};  nW = 1e5;
[Aq, bq] = buildSystem('ghost3', uq);
for r = [1 5 9]
    i = mod(r-1,3)+1;  j = (r-i)/3 + 1;
    [est, se] = walkGrid(Aq, bq, r, nW);
    fprintf('  grid walk (%d,%d): %.3f +- %.3f (exact %g)\n', i, j, est, 1.96*se, uq(i, j));
    nFail = nFail + ~check(sprintf('grid walk (%d,%d) within 4 standard errors', i, j), abs(est - uq(i, j)) < 4*se);
    [est, se] = walkOnSpheres(i, j, nW, 1e-4, uq);
    fprintf('  walk-on-spheres (%d,%d): %.3f +- %.3f (exact %g)\n', i, j, est, 1.96*se, uq(i, j));
    nFail = nFail + ~check(sprintf('walk-on-spheres (%d,%d) within 4 standard errors', i, j), abs(est - uq(i, j)) < 4*se + 1e-3);
end

% ---- 8. physics / interpolation sanity
disp('Test 8: maximum principle and Lagrange weights');
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
