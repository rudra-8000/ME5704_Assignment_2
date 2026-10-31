% A5_MONTECARLO  Random-walk solutions: (i) walk on the h=1 cubic-ghost (M3) system,
%   (ii) mesh-free walk-on-spheres.  Both are compared with deterministic references.
%   (i) must reproduce the M3 answer (same equations);  (ii) estimates the true
%   continuum temperature, i.e. the discretisation-free benchmark.
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
[A, b] = buildSystem('ghost3');  T3 = gaussElim(A, b);
load(fullfile(resDir, 'a2_refinement.mat'), 'Tref');
Tref = reshape(Tref, [], 1);                          % same node ordering r = (j-1)*3 + i
epsAbs = 1e-4;

% ---- error vs number of walks
Ns = round(logspace(2, 5, 7));  rmsG = zeros(size(Ns));  rmsW = zeros(size(Ns));  reps = 5;
for k = 1:numel(Ns)
    eG = zeros(9, reps);  eW = zeros(9, reps);        % errors of 9 nodes x 5 independent repeats
    for r = 1:9
        i = mod(r-1, 3) + 1;  j = (r - i)/3 + 1;
        for q = 1:reps
            eG(r,q) = walkGrid(A, b, r, Ns(k)) - T3(r);
            eW(r,q) = walkOnSpheres(i, j, Ns(k), epsAbs) - Tref(r);
        end
    end
    rmsG(k) = sqrt(mean(eG(:).^2));  rmsW(k) = sqrt(mean(eW(:).^2));
    fprintf('N = %6d   rms error  grid walk %.3f K   walk-on-spheres %.3f K\n', Ns(k), rmsG(k), rmsW(k));
end

% ---- final table with 95% confidence intervals
N = 2e5;  est = zeros(9,2);  se = zeros(9,2);
for r = 1:9
    i = mod(r-1, 3) + 1;  j = (r - i)/3 + 1;
    [est(r,1), se(r,1)] = walkGrid(A, b, r, N);
    [est(r,2), se(r,2)] = walkOnSpheres(i, j, N, epsAbs);
end
inCI = abs(est(:,2) - Tref) < 1.96*se(:,2);
fprintf('\n N = %d walks per node\n node   M3      gridWalk  +-95%%   Tref    WoS    +-95%%\n', N);
for r = 1:9
    fprintf('(%d,%d) %7.2f %8.2f %6.2f %8.2f %7.2f %6.2f\n', mod(r-1,3)+1, ceil(r/3), T3(r), est(r,1), 1.96*se(r,1), Tref(r), est(r,2), 1.96*se(r,2));
end
fprintf('WoS reference inside the 95%% CI at %d of 9 nodes; max |WoS - Tref| = %.2f K\n', sum(inCI), max(abs(est(:,2) - Tref)));
save(fullfile(resDir, 'a5_montecarlo.mat'), 'Ns', 'rmsG', 'rmsW', 'est', 'se', 'N');

% ---- Figure: error decay
f = figure('Visible','off');
loglog(Ns, rmsG, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);  hold on;
loglog(Ns, rmsW, 's-', 'LineWidth', 1.2, 'MarkerSize', 4);
loglog(Ns, 60./sqrt(Ns), 'k--');  text(Ns(2), 90/sqrt(Ns(2)), '$N^{-1/2}$', 'Interpreter', 'latex');
xlabel('number of walks $N$', 'Interpreter', 'latex');  ylabel('rms error [K]');  grid on;  set(gca, 'FontSize', 9);
legend('grid walk vs M3', 'walk-on-spheres vs $T^*$', 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 7);
saveFig(f, 'a_montecarlo', 8, 6.5);
