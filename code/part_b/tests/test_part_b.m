% TEST_PART_B  Validation of the part (b) scripts against known answers.
%   1. Lecture Notes Sec. 5.5.2 example: Richardson on f = 3^x/x at x = 1       -> fdWeights
%   2. Exactness on the h = 1 grid: quadratic (O(h^2) schemes exact), cubic (C4 exact) -> gradGiven
%   3. Order of accuracy on f = e^x sin y, h halved four times                    -> fdWeights, fdNonUniform
%   4. FD-step check on the finest refined field                                   -> gradFine
%   5. Edge heat rate on functions with a known boundary integral                  -> edgeFlux
%   6. Convergence of Q/k with the grid                                            -> edgeFlux on refined fields
%   7. MATLAB's gradient (used ONLY here, as an independent oracle)                -> gradFine
nFail = 0;
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
P = plateProblem();

% ---- 1. lecture example (notes: D(0.2) = 0.3298, D(0.1) = 0.3078, D = 0.3005, exact 3 ln 3 - 3 = 0.2958)
disp('Test 1: lecture example, Richardson on f = 3^x/x at x = 1');
f = @(x) 3.^x ./ x;
Dh = @(h) fdApply(f, 1, 'F2', h);
D1 = Dh(0.2);  D2 = Dh(0.1);  Dr = 4/3*D2 - 1/3*D1;
fprintf('  D(0.2) = %.4f  D(0.1) = %.4f  Richardson %.4f  exact %.4f\n', D1, D2, Dr, 3*log(3) - 3);
nFail = nFail + ~check('matches the notes to 4 decimals', ...
    all(abs(round([D1 D2 Dr], 4) - [0.3298 0.3078 0.3005]) < 1e-12));

% ---- 2. exactness on the given grid (edge and line values taken from u)
disp('Test 2: exactness on the h = 1 grid');
[Xn, Yn] = ndgrid(1:3, 1:3);
uq = @(x,y) x.^2 - y.^2 + x.*y;   uqx = @(x,y) 2*x + y;      uqy = @(x,y) -2*y + x;
uc = @(x,y) x.^3 - 3*x.*y.^2;     ucx = @(x,y) 3*x.^2 - 3*y.^2;  ucy = @(x,y) -6*x.*y;
Dq = gradGiven(uq(Xn, Yn), uq);   Dc = gradGiven(uc(Xn, Yn), uc);
eq = @(D, s, gx, gy) max(abs([D.(s).Tx(:) - gx(Xn(:),Yn(:)); D.(s).Ty(:) - gy(Xn(:),Yn(:))]), [], 'omitnan');
fprintf('  quadratic: C2 %.1e  OS2 %.1e  NU %.1e  C4 %.1e  R %.1e  | OS1 %.2f (O(h), not exact)\n', ...
    eq(Dq,'C2',uqx,uqy), eq(Dq,'OS2',uqx,uqy), eq(Dq,'NU',uqx,uqy), eq(Dq,'C4',uqx,uqy), eq(Dq,'R',uqx,uqy), eq(Dq,'OS1',uqx,uqy));
fprintf('  cubic:     C4 %.1e  R %.1e  | C2 %.2f (O(h^2), not exact)\n', eq(Dc,'C4',ucx,ucy), eq(Dc,'R',ucx,ucy), eq(Dc,'C2',ucx,ucy));
nFail = nFail + ~check('quadratic: C2, OS2, NU, C4, R exact (< 1e-10)', ...
    max([eq(Dq,'C2',uqx,uqy), eq(Dq,'OS2',uqx,uqy), eq(Dq,'NU',uqx,uqy), eq(Dq,'C4',uqx,uqy), eq(Dq,'R',uqx,uqy)]) < 1e-10);
nFail = nFail + ~check('cubic: C4 and R exact, C2 not', ...
    max(eq(Dc,'C4',ucx,ucy), eq(Dc,'R',ucx,ucy)) < 1e-10 && eq(Dc,'C2',ucx,ucy) > 0.1);

% ---- 3. order of accuracy on f = e^x sin y, derivative in x at (1,1)
disp('Test 3: order of accuracy, f = e^x sin y');
g = @(x) exp(x) * sin(1);  exact = exp(1) * sin(1);
hs = 0.4 ./ 2.^(0:4);
names = {'F1', 'B1', 'C2', 'F2', 'B2', 'NU', 'C4', 'R'};  expected = [1 1 2 2 2 2 4 4];
fid = fopen(fullfile(resDir, 'test_order.csv'), 'w');
fprintf(fid, 'scheme,expected,err_first,err_last,observed\n');
for q = 1:numel(names)
    e = zeros(size(hs));
    for k = 1:numel(hs)
        h = hs(k);
        switch names{q}
            case 'NU', d = fdNonUniform(h, h/2) * [g(1-h); g(1); g(1+h/2)];      % spacing h left, h/2 right
            case 'R',  d = 4/3*fdApply(g, 1, 'C2', h) - 1/3*fdApply(g, 1, 'C2', 2*h);
            otherwise, d = fdApply(g, 1, names{q}, h);
        end
        e(k) = abs(d - exact);
    end
    obs = log2(e(end-1) / e(end));
    fprintf('  %-3s errors %s  observed order %.2f (expected %d)\n', names{q}, mat2str(e, 3), obs, expected(q));
    fprintf(fid, '%s,%d,%.6e,%.6e,%.4f\n', names{q}, expected(q), e(1), e(end), obs);
    nFail = nFail + ~check(sprintf('%s has order %d (+-0.15)', names{q}, expected(q)), abs(obs - expected(q)) < 0.15);
end
fclose(fid);

% ---- 4-7 use the refined fields of b2_reference
B2 = load(fullfile(resDir, 'b2_reference.mat'));
disp('Test 4: FD-step check on the finest field');
fprintf('  Richardson in the step changes grad T by %.4f K/m\n', B2.stepChange);
nFail = nFail + ~check('step error on the N = 256 field below 0.05 K/m', B2.stepChange < 0.05);

disp('Test 5: edge heat rate on known boundary integrals');
for hq = [1 0.5]
    n = P.L/hq;  [X, Y] = ndgrid((0:n)*hq);
    Q1 = edgeFlux(X.^2, hq, 2, false);                % -(boundary integral of dT/dn) = -(8 * 4) = -32
    Q2 = edgeFlux(X.^2 - Y.^2 + X.*Y, hq, 2, false);  % zero (Laplacian is zero)
    fprintf('  h = %.1f:  T = x^2 -> %.12f (exact -32);  T = x^2-y^2+xy -> %.1e (exact 0)\n', hq, Q1, Q2);
    nFail = nFail + ~check(sprintf('edgeFlux exact for quadratics at h = %.1f', hq), abs(Q1 + 32) < 1e-10 && abs(Q2) < 1e-10);
end

disp('Test 6: convergence of Q/k with the grid');
dQ = -diff(B2.Qk);  ratios = dQ(1:end-1) ./ dQ(2:end);
fprintf('  Q/k = %s, ratio of successive changes %s\n', mat2str(B2.Qk, 6), mat2str(ratios, 3));
nFail = nFail + ~check('Q/k decreases monotonically, first order (ratio 1.6 to 2.4)', ...
    all(dQ > 0) && all(ratios > 1.6 & ratios < 2.4));

disp('Test 7: gradFine vs MATLAB gradient (oracle)');
n = numel(B2.Ns);  hn = P.L/B2.Ns(n);
[FY, FX] = gradient(B2.Tfine, hn);                    % T(i,j): i runs along x (dim 1) -> second output
[gx, gy] = gradFine(B2.Tfine, B2.Ns(n), 1);
idx = (1:3) * B2.Ns(n)/P.L + 1;
dev = max(abs([gx(:) - reshape(FX(idx,idx), [], 1); gy(:) - reshape(FY(idx,idx), [], 1)]));
fprintf('  max |gradFine - gradient| = %.1e\n', dev);
nFail = nFail + ~check('agrees to round-off', dev < 1e-9);

fprintf('\n%d test(s) failed\n', nFail);
if nFail > 0, error('part (b) tests failed'); end

function d = fdApply(fun, x, name, h)
% Apply the Table 3-5 scheme NAME to the function FUN at X with step H.
[off, w] = fdWeights(name, h);
d = w * reshape(fun(x + off*h), [], 1);
end

function ok = check(name, cond)
% Print one PASS/FAIL line and return the result.
label = 'FAIL';  if cond, label = 'PASS'; end
fprintf('  [%s] %s\n', label, name);  ok = cond;
end
