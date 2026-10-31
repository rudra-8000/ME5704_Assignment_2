% A6_ROBUSTNESS  Does the ranking of the h = 1 treatments depend on where the source line ends?
%   The h = 1 systems only see that the line passes the node rows y = 2 and 3, so they are the
%   same for every line below.  The converged field is not.  For each line the reference T* is
%   rebuilt (grids N = 32, 64, 128, Richardson with the observed order, as in a2_refinement)
%   and the three line-based treatments (M1 to M3) are compared with it.
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
lines = [2 3; 1.5 3.5; 1.5 3; 2 3.5];                % [yLo yHi]; the first is the given problem
meth  = {'sw', 'ghost1', 'ghost3'};
Ns = [32 64 128];  dev = zeros(size(lines,1), 3);  pMean = zeros(size(lines,1), 1);
for g = 1:size(lines,1)
    S = zeros(numel(Ns), 3, 3);
    for k = 1:numel(Ns)
        T = solveFine(Ns(k), [], lines(g,1), lines(g,2));
        idx = (1:3) * Ns(k)/4 + 1;  S(k,:,:) = T(idx, idx);
    end
    d1 = squeeze(S(1,:,:) - S(2,:,:));  d2 = squeeze(S(2,:,:) - S(3,:,:));
    p = log2(d1 ./ d2);  Tstar = squeeze(S(3,:,:)) - d2 ./ (2.^p - 1);
    pMean(g) = mean(p(:));
    for m = 1:3
        [A, b] = buildSystem(meth{m});
        dev(g,m) = max(max(abs(reshape(gaussElim(A, b), 3, 3) - Tstar)));
    end
    fprintf('line y in [%g, %g]: mean observed p %.2f;  max dev. M1 %.1f, M2 %.1f, M3 %.1f K\n', lines(g,:), pMean(g), dev(g,:));
end
save(fullfile(resDir, 'a6_robustness.mat'), 'lines', 'dev', 'pMean');
