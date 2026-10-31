% A3_SOLVERS  Same linear system A*T = b (cubic ghost node M3, h = 1), different solvers.
%   Direct: Gauss elimination, Crout LU.  Iterative: Jacobi, Gauss-Seidel (Liebmann), SOR.
[A, b] = buildSystem('ghost3');
n = numel(b);  tol = 1e-8;
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');

% ---- diagonal dominance (sufficient condition for Jacobi / Gauss-Seidel convergence)
offSum = sum(abs(A), 2) - abs(diag(A));
fprintf('diagonally dominant: %d (strict rows: %d of %d)\n', all(abs(diag(A)) >= offSum), sum(abs(diag(A)) > offSum), n);

% ---- direct solvers
xG = gaussElim(A, b);  xL = luSolve(A, b);
fprintf('max |Gauss - LU| = %.2e\n', max(abs(xG - xL)));

% ---- iterative solvers: sweeps to reach residual < tol
[xJ, hJ]  = jacobi(A, b, tol, 5000);
[xGS, hGS] = sor(A, b, 1.0, tol, 5000);
omegas = 1.0:0.05:1.95;  nit = zeros(size(omegas));
for k = 1:numel(omegas)
    [~, hh] = sor(A, b, omegas(k), tol, 5000);  nit(k) = numel(hh);
end
[~, kBest] = min(nit);  omegaOpt = omegas(kBest);
[xS, hS] = sor(A, b, omegaOpt, tol, 5000);
fprintf('sweeps to residual < %.0e:  Jacobi %d, Gauss-Seidel %d, SOR(omega=%.2f) %d\n', tol, numel(hJ), numel(hGS), omegaOpt, numel(hS));
rate = @(h) (h(numel(h))/h(ceil(numel(h)/2)))^(1/(numel(h) - ceil(numel(h)/2)));   % mean residual factor per sweep (2nd half)
rates = [rate(hJ) rate(hGS) rate(hS)];  sweeps = [numel(hJ) numel(hGS) numel(hS)];
fprintf('mean residual factor per sweep: Jacobi %.3f, GS %.3f, SOR %.3f  (Jacobi^2 = %.3f)\n', rate(hJ), rate(hGS), rate(hS), rate(hJ)^2);
fprintf('max |iterative - Gauss|: Jacobi %.1e, GS %.1e, SOR %.1e\n', max(abs(xJ-xG)), max(abs(xGS-xG)), max(abs(xS-xG)));

% ---- how the sweep count scales with the grid for GS vs SOR (matrix-free, fine grids)
Ns = [8 16 32];  gsSweeps = zeros(size(Ns));  sorSweeps = zeros(size(Ns));
for k = 1:numel(Ns)
    [~, sorSweeps(k)] = solveFine(Ns(k));
    [~, gsSweeps(k)] = solveFine(Ns(k), 1.0);
end
fprintf('fine-grid sweeps  N = %s:  GS %s | SOR(opt) %s\n', mat2str(Ns), mat2str(gsSweeps), mat2str(sorSweeps));
save(fullfile(resDir, 'a3_solvers.mat'), 'hJ', 'hGS', 'hS', 'omegas', 'nit', 'omegaOpt', 'Ns', 'gsSweeps', 'sorSweeps', 'rates', 'sweeps');

% ---- Figure: residual histories + sweeps vs omega
f = figure('Visible','off');
subplot(1,2,1);  semilogy(1:numel(hJ), hJ, 'LineWidth', 1.2);  hold on;
semilogy(1:numel(hGS), hGS, 'LineWidth', 1.2);  semilogy(1:numel(hS), hS, 'LineWidth', 1.2);
xlabel('sweep');  ylabel('max residual');  grid on;  set(gca, 'FontSize', 8);
legend('Jacobi', 'Gauss-Seidel', sprintf('SOR %.2f', omegaOpt), 'Location', 'northeast', 'FontSize', 7);
subplot(1,2,2);  plot(omegas, nit, 'o-', 'LineWidth', 1.2, 'MarkerSize', 3);
xlabel('relaxation factor $\omega$', 'Interpreter', 'latex');  ylabel('sweeps to $10^{-8}$', 'Interpreter', 'latex');  grid on;  set(gca, 'FontSize', 8);
saveFig(f, 'a_solvers', 15, 5.5);
