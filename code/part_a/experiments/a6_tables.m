% A6_TABLES  Turn the saved results into complete LaTeX tabulars (report/tables/*.tex), so that
%   no number in the report is typed by hand.  The report \input{}s each file inside a table float.
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
tabDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'report', 'tables');
if ~exist(tabDir, 'dir'), mkdir(tabDir); end
R  = load(fullfile(resDir, 'a1_treatments.mat'));
R2 = load(fullfile(resDir, 'a2_refinement.mat'));
R3 = load(fullfile(resDir, 'a3_solvers.mat'));
R4 = load(fullfile(resDir, 'a4_transient.mat'));
R5 = load(fullfile(resDir, 'a5_montecarlo.mat'));
P  = plateProblem();
meth = {'sw','ghost1','ghost3','snapL','snapR','snapAll'};

% ---- Table 1: nine node temperatures for every treatment + reference
fid = openTab(tabDir, 'a_treatments', '@{}l rrrrrr r@{}', ...
    'Node & M1 & M2 & M3 & M4 & M5 & M6 & $T^*$');
for j = 3:-1:1, for i = 1:3
    fprintf(fid, '$(%d,%d)$', i, j);
    for m = 1:6, fprintf(fid, ' & %.1f', R.(meth{m})(i,j)); end
    fprintf(fid, ' & %.1f \\\\\n', R2.Tref(i,j));
end, end
fprintf(fid, '\\midrule\nmax dev. [K]');
for m = 1:6, fprintf(fid, ' & %.1f', max(max(abs(R.(meth{m}) - R2.Tref)))); end
fprintf(fid, ' & n/a \\\\\nmax dev. [\\%%]');
for m = 1:6, fprintf(fid, ' & %.1f', 100*max(max(abs(R.(meth{m}) - R2.Tref) ./ R2.Tref))); end
fprintf(fid, ' & n/a \\\\\n');  closeTab(fid);

% ---- Answer table: M3 temperatures laid out like the grid, with T* in brackets
fid = openTab(tabDir, 'a_answer', '@{}c|ccc@{}', '$y \backslash x$ & $x=1$ & $x=2$ & $x=3$');
for j = 3:-1:1
    fprintf(fid, '$y=%d$', j);
    for i = 1:3, fprintf(fid, ' & \\textbf{%.1f} [%.1f]', R.ghost3(i,j), R2.Tref(i,j)); end
    fprintf(fid, ' \\\\\n');
end
closeTab(fid);

% ---- Table 2: refinement sequence at three representative nodes
fid = openTab(tabDir, 'a_refinement', '@{}r r rrr r@{}', '$N$ & $h$ & $T_{2,2}$ & $T_{1,3}$ & $T_{3,1}$ & SOR sweeps');
for k = 1:numel(R2.Ns)
    fprintf(fid, '$%d$ & %g & %.2f & %.2f & %.2f & %d \\\\\n', R2.Ns(k), P.L/R2.Ns(k), ...
            R2.S(k,2,2), R2.S(k,1,3), R2.S(k,3,1), R2.iters(k));
end
fprintf(fid, '\\midrule\nobserved $p$ & & %.2f & %.2f & %.2f & \\\\\n', R2.pObs(2,2), R2.pObs(1,3), R2.pObs(3,1));
fprintf(fid, 'Richardson $T^*$ & & %.2f & %.2f & %.2f & \\\\\n', R2.Tref(2,2), R2.Tref(1,3), R2.Tref(3,1));
closeTab(fid);

% ---- Table 3: solvers on the 9x9 system
fid = openTab(tabDir, 'a_solvers', '@{}l l r r@{}', 'Method & Type & Sweeps to $10^{-8}$ & Residual factor per sweep');
fprintf(fid, 'Gauss elimination (pivoting) & direct & n/a & n/a \\\\\nCrout LU & direct & n/a & n/a \\\\\n');
nm = {'Jacobi', 'Gauss--Seidel', sprintf('SOR, $\\omega=%.2f$', R3.omegaOpt)};
for q = 1:3, fprintf(fid, '%s & iterative & %d & %.3f \\\\\n', nm{q}, R3.sweeps(q), R3.rates(q)); end
closeTab(fid);

% ---- Table 4: time marching schemes
[A, b] = buildSystem('ghost3');  T3 = gaussElim(A, b);                  % steady value of T_{2,2} for the % column
fid = openTab(tabDir, 'a_transient', '@{}l r r r r@{}', 'Scheme & $\Delta t/\Delta t_{\max}$ & Steps & $\max|\Delta T_{2,2}|$ [K] & [\%]');
nm = struct('euler', 'Explicit Euler', 'rk4', 'RK4', 'ieuler', 'Implicit Euler', 'cn', 'Crank--Nicolson');
for k = 1:size(R4.runs, 1)
    fprintf(fid, '%s & %.1f & %d & %.2f & %.2f \\\\\n', nm.(R4.runs{k,1}), R4.tab(k,1), R4.tab(k,2), R4.tab(k,3), 100*R4.tab(k,3)/T3(5));
end
closeTab(fid);

% ---- Table 5: Monte Carlo at N walks per node, with 95% half-widths
Tref = reshape(R2.Tref, [], 1);
fid = openTab(tabDir, 'a_montecarlo', '@{}l r r r r@{}', 'Node & M3 & Grid walk $\pm$95\% & $T^*$ & Walk-on-spheres $\pm$95\%');
for r = 9:-1:1
    i = mod(r-1, 3) + 1;  j = (r - i)/3 + 1;
    fprintf(fid, '$(%d,%d)$ & %.2f & %.2f $\\pm$ %.2f & %.2f & %.2f $\\pm$ %.2f \\\\\n', i, j, T3(r), R5.est(r,1), 1.96*R5.se(r,1), Tref(r), R5.est(r,2), 1.96*R5.se(r,2));
end
closeTab(fid);

% ---- Table 6: convergence test (exact harmonic function)
C = readCsv(fullfile(resDir, 'test_convergence.csv'));
fid = openTab(tabDir, 'a_test', '@{}c r r@{}', '$h$ & $\max|T-u|$ & observed order');
for k = 1:size(C,1)
    ord = 'n/a';  if k > 1, ord = sprintf('%.2f', log2(C(k-1,2)/C(k,2))); end
    fprintf(fid, '$1/%d$ & %.2e & %s \\\\\n', round(1/C(k,1)), C(k,2), ord);
end
closeTab(fid);

function fid = openTab(tabDir, name, spec, header)
% Open NAME.tex and write the tabular preamble (column spec, top rule, header row, mid rule).
fid = fopen(fullfile(tabDir, [name '.tex']), 'w');
fprintf(fid, '\\begin{tabular}{%s}\n\\toprule\n%s \\\\\n\\midrule\n', spec, header);
end

function closeTab(fid)
% Write the bottom rule and close the tabular.
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');  fclose(fid);
end

function M = readCsv(file)
% Read a two-column numeric CSV with one header line.
fid = fopen(file);  fgetl(fid);  M = fscanf(fid, '%f,%f', [2 Inf])';  fclose(fid);
end
