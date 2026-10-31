% A2_REFINEMENT  Grid refinement (line on nodes) + Richardson extrapolation.
%   Solves the plate on N = 8..128 divisions, samples the 9 given nodes, estimates the
%   observed order p from the three finest grids and extrapolates to a reference field.
%   Then measures the error of every h = 1 treatment against that reference.
P = plateProblem();
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
Ns = [8 16 32 64 128];  S = zeros(numel(Ns), 3, 3);  iters = zeros(size(Ns));
for k = 1:numel(Ns)
    [T, iters(k)] = solveFine(Ns(k));
    idx = (1:3) * Ns(k)/4 + 1;                       % indices of x,y = 1,2,3
    S(k,:,:) = T(idx, idx);
    if Ns(k) == 128, Tfine = T; end                  % keep for the contour plot
end
% ---- observed order  p = log2( (T_h - T_h/2) / (T_h/2 - T_h/4) )  from the 3 finest grids
d1 = squeeze(S(end-2,:,:) - S(end-1,:,:));  d2 = squeeze(S(end-1,:,:) - S(end,:,:));
pObs = log2(d1 ./ d2);
% ---- Richardson:  T* = T_h/4 + (T_h/4 - T_h/2) / (2^p - 1),  node-wise p
Tfin  = squeeze(S(end,:,:));  Tprev = squeeze(S(end-1,:,:));
Tref  = Tfin + (Tfin - Tprev) ./ (2.^pObs - 1);
Tref2 = Tfin + (Tfin - Tprev) / 3;                   % textbook assumption p = 2
fprintf('observed order p (rows x=1..3, cols y=1..3):\n'); disp(pObs);
fprintf('Richardson reference (p observed), rows y=3..1:\n'); disp(flipud(Tref'));
fprintf('Richardson with p=2 assumed, rows y=3..1:\n');       disp(flipud(Tref2'));
fprintf('SOR sweeps per grid: '); fprintf('%d ', iters); fprintf('\n');

% ---- check of the reference: repeat the extrapolation one level finer (N = 64, 128, 256)
T256 = solveFine(256);  idx = (1:3) * 64 + 1;  S256 = T256(idx, idx);
S64 = squeeze(S(end-1,:,:));  S128 = squeeze(S(end,:,:));
pFine = log2((S64 - S128) ./ (S128 - S256));
TrefFine = S256 + (S256 - S128) ./ (2.^pFine - 1);
dev256 = max(max(abs(TrefFine - Tref)));
fprintf('reference from N = 64,128,256 differs from the one from N = 32,64,128 by at most %.3f K (N = 256 itself: %.3f K)\n', dev256, max(max(abs(S256 - Tref))));

% ---- errors of the h=1 treatments against the reference
R = load(fullfile(resDir, 'a1_treatments.mat'));
fn = fieldnames(R);  errMax = zeros(numel(fn), 1);
for m = 1:numel(fn)
    errMax(m) = max(max(abs(R.(fn{m}) - Tref)));
    fprintf('%-8s  max|T - Tref| = %6.2f K\n', fn{m}, errMax(m));
end
save(fullfile(resDir, 'a2_refinement.mat'), 'Ns', 'S', 'pObs', 'Tref', 'Tref2', 'Tfine', 'iters', 'fn', 'errMax', 'dev256');
fid = fopen(fullfile(resDir, 'part_a_temperature.csv'), 'w');
fprintf(fid, 'i,j,x,y,T_cubic_h1,T_sw_h1,T_reference\n');   % interface for parts (b), (c)
for j = 1:3, for i = 1:3
    fprintf(fid, '%d,%d,%d,%d,%.4f,%.4f,%.4f\n', i, j, i, j, R.ghost3(i,j), R.sw(i,j), Tref(i,j));
end, end
fclose(fid);

% ---- Figure: grid convergence from successive differences |T_h - T_h/2| (no reference needed)
f = figure('Visible','off');  hold on;  h = P.L ./ Ns;
D = abs(S(1:end-1,:,:) - S(2:end,:,:));  D = reshape(D, numel(Ns)-1, 9);
loglog(h(1:end-1), D(:,1:9), '-o', 'Color', [0.55 0.55 0.55], 'MarkerSize', 3, 'HandleVisibility', 'off');
hn = loglog(h(1:end-1), D(:,5), 'o-', 'Color', [0 0.45 0.74], 'LineWidth', 1.4, 'MarkerSize', 4);   % node (2,2)
h1 = loglog(h(1:end-1), D(1,5)*h(1:end-1)/h(1), 'k--');  h2 = loglog(h(1:end-1), D(1,5)*(h(1:end-1)/h(1)).^2, 'k:');
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 9);
xlabel('grid spacing $h$', 'Interpreter', 'latex');  ylabel('$|T_h - T_{h/2}|$ [K]', 'Interpreter', 'latex');  grid on;  box on;
legend([hn h1 h2], {'node (2,2); grey: other eight', '$O(h)$', '$O(h^2)$'}, 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 7);
saveFig(f, 'a_convergence', 8, 6.5);

% ---- Figure: converged temperature field (N = 128)
f = figure('Visible','off');  x = linspace(0, P.L, 129);
imagesc(x, x, Tfine', [200 350]);  set(gca, 'YDir', 'normal');  hold on;  colormap(turbo);
cb = colorbar;  cb.Label.String = 'T [K]';
plot([P.xLine P.xLine], [P.yLo P.yHi], 'k-', 'LineWidth', 2);
plot(repmat(1:3, 1, 3), repelem(1:3, 3), 'wo', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
axis equal tight;  xlabel('x');  ylabel('y');  set(gca, 'FontSize', 9);
saveFig(f, 'a_field', 8, 6.5);
