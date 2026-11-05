% B2_REFERENCE  Reference gradient and heat rate from part (a)'s refined fields (line on nodes).
%   For N = 32..256: field from solveFine(N), central gradient at the nine node positions
%   (gradFine) and the edge heat rate Q/k (edgeFlux). Then, as for T* in part (a):
%   observed order p from the three finest grids and Richardson extrapolation.
P = plateProblem();
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
Ns = [32 64 128 256];  nN = numel(Ns);
GX = zeros(nN, 3, 3);  GY = zeros(nN, 3, 3);  Qk = zeros(1, nN);  QE = zeros(nN, 4);
for k = 1:nN
    T = solveFine(Ns(k));
    [gx, gy] = gradFine(T, Ns(k), 1);
    GX(k,:,:) = gx;  GY(k,:,:) = gy;
    [Qk(k), QE(k,:)] = edgeFlux(T, P.L/Ns(k));
    fprintf('N = %3d  h = %.4f  Q/k = %.3f K m\n', Ns(k), P.L/Ns(k), Qk(k));
    if k == nN, Tfine = T;  [gx2, gy2] = gradFine(T, Ns(k), 2); end
end

% ---- FD-step check on the finest field: Richardson in the step (Sec. 5.5.2), D = 4/3 D(h) - 1/3 D(2h)
stepChange = max(abs([4/3*gx - 1/3*gx2 - gx, 4/3*gy - 1/3*gy2 - gy]), [], 'all');
fprintf('N = %d: Richardson in the FD step changes grad T by at most %.4f K/m\n', Ns(end), stepChange);

% ---- grid refinement: observed order and Richardson (same procedure as T* in part a)
[GXref, pX] = extrapolate(GX(end-2:end,:,:));
[GYref, pY] = extrapolate(GY(end-2:end,:,:));
[Qref,  pQ] = extrapolate(Qk(end-2:end)');
QEref = zeros(1, 4);  for e = 1:4, QEref(e) = extrapolate(QE(end-2:end, e)); end
% same extrapolation one level coarser (N = 32, 64, 128): spread = uncertainty of the reference
GXc = extrapolate(GX(1:3,:,:));  GYc = extrapolate(GY(1:3,:,:));  Qc = extrapolate(Qk(1:3)');
refSpread = max(abs([GXref(:) - GXc(:); GYref(:) - GYc(:)]));
fprintf('observed order p, dT/dx (rows x = 1..3, cols y = 1..3):\n');  disp(pX);
fprintf('observed order p, dT/dy:\n');  disp(pY);
fprintf('reference dT/dx, rows y = 3..1:\n');  disp(flipud(GXref'));
fprintf('reference dT/dy, rows y = 3..1:\n');  disp(flipud(GYref'));
fprintf('reference from N = 32,64,128 differs from N = 64,128,256 by at most %.2f K/m\n', refSpread);
fprintf('Q/k: observed p = %.2f, Richardson Q*/k = %.2f K m (from N = 32,64,128: %.2f)\n', pQ, Qref, Qc);
save(fullfile(resDir, 'b2_reference.mat'), 'Ns', 'GX', 'GY', 'Qk', 'QE', 'GXref', 'GYref', 'pX', 'pY', ...
     'Qref', 'pQ', 'QEref', 'Qc', 'refSpread', 'stepChange', 'Tfine');

% ---- Figure: (left) change of grad T between successive grids, (right) Q/k against h
f = figure('Visible', 'off');  tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
h = P.L ./ Ns;  hh = h(1:end-1);
nexttile;  hold on;
D = abs([reshape(GX(1:end-1,:,:) - GX(2:end,:,:), nN-1, 9), reshape(GY(1:end-1,:,:) - GY(2:end,:,:), nN-1, 9)]);
loglog(hh, D, '-o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3, 'HandleVisibility', 'off');
h1 = loglog(hh, max(D(1,:))*hh/hh(1), 'k--');  h2 = loglog(hh, max(D(1,:))*(hh/hh(1)).^2, 'k:');
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 8);  grid on;  box on;
xlabel('grid spacing $h$ [m]', 'Interpreter', 'latex');  ylabel('$|\nabla T_h-\nabla T_{h/2}|$ [K/m]', 'Interpreter', 'latex');
legend([h1 h2], {'$O(h)$', '$O(h^2)$'}, 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 7);
title('18 gradient components', 'FontWeight', 'normal', 'FontSize', 8);
nexttile;  hold on;
plot(h, Qk, 'o-', 'Color', [0 0.45 0.74], 'LineWidth', 1.2, 'MarkerSize', 4);
yline(Qref, 'k--', sprintf('Richardson %.1f', Qref), 'FontSize', 7, 'LabelHorizontalAlignment', 'right');
set(gca, 'FontSize', 8);  grid on;  box on;  xlim([0 0.14]);
xlabel('grid spacing $h$ [m]', 'Interpreter', 'latex');  ylabel('$Q/k$ [K\,m]', 'Interpreter', 'latex');
title('heat rate through the edges', 'FontWeight', 'normal', 'FontSize', 8);
saveFigB(f, 'b_convergence', 16, 6.5);

function [ref, p] = extrapolate(S)
% S(1..3, ...) = values on grids h, h/2, h/4.  p = log2((S1-S2)/(S2-S3)),  ref = S3 + (S3-S2)/(2^p-1).
d1 = S(1,:) - S(2,:);  d2 = S(2,:) - S(3,:);
assert(all(d1 ./ d2 > 0), 'extrapolate: refinement sequence is not monotone');
p   = log2(d1 ./ d2);
ref = S(3,:) + (S(3,:) - S(2,:)) ./ (2.^p - 1);
sz  = size(S);  if numel(sz) > 2, p = reshape(p, sz(2:end));  ref = reshape(ref, sz(2:end)); end
end
