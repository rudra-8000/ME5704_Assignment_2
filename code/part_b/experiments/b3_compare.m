% B3_COMPARE  Error of every h = 1 scheme against the reference, split into
%     data error            = scheme on T_cubic_h1  minus  same scheme on T_reference   (part a's error)
%     differentiation error = scheme on T_reference minus  reference gradient           (h = 1 spacing)
%   and the choice of the headline scheme: in each group of nodes that allows the same schemes,
%   the candidate with the smallest largest error against the reference.
P = plateProblem();
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
B1 = load(fullfile(resDir, 'b1_given_grid.mat'));
B2 = load(fullfile(resDir, 'b2_reference.mat'));
ref = struct('Tx', B2.GXref, 'Ty', B2.GYref);

% groups: component, nodes [i j], candidate schemes (C2 at the line nodes is the baseline N0)
G = struct('comp', {'Tx', 'Tx', 'Tx', 'Tx', 'Ty', 'Ty'}, ...
           'nodes', {[1 2; 1 3], [2 2; 2 3], [2 1], [1 1; 3 1; 3 2; 3 3], [1 2; 2 2; 3 2], [1 1; 2 1; 3 1; 1 3; 2 3; 3 3]}, ...
           'cand', {{'C2','OS1','NU'}, {'C2','OS1','OS2','NU'}, {'C2','C4'}, {'C2'}, {'C2','C4'}, {'C2'}});
Tx = nan(3);  Ty = nan(3);  schX = cell(3);  schY = cell(3);
for g = 1:numel(G)
    idx = sub2ind([3 3], G(g).nodes(:,1), G(g).nodes(:,2));
    r = ref.(G(g).comp)(idx);
    nc = numel(G(g).cand);  G(g).err = zeros(1,nc);  G(g).dataErr = zeros(1,nc);  G(g).diffErr = zeros(1,nc);
    for c = 1:nc
        s  = G(g).cand{c};
        vC = B1.Dcub.(s).(G(g).comp)(idx);  vR = B1.Dref.(s).(G(g).comp)(idx);
        G(g).err(c)     = max(abs(vC - r));
        G(g).dataErr(c) = max(abs(vC - vR));
        G(g).diffErr(c) = max(abs(vR - r));
    end
    [~, best] = min(G(g).err);  G(g).best = G(g).cand{best};
    v = B1.Dcub.(G(g).best).(G(g).comp)(idx);
    if strcmp(G(g).comp, 'Tx'), Tx(idx) = v;  schX(idx) = {G(g).best};
    else,                      Ty(idx) = v;  schY(idx) = {G(g).best}; end
    fprintf('group %d (%s at %s): ', g, G(g).comp, mat2str(G(g).nodes));
    for c = 1:nc, fprintf('%s %.2f (data %.2f, diff %.2f)  ', G(g).cand{c}, G(g).err(c), G(g).dataErr(c), G(g).diffErr(c)); end
    fprintf('-> %s\n', G(g).best);
end
assert(~any(isnan([Tx(:); Ty(:)])), 'b3_compare: a node has no headline value');

qmag = hypot(Tx, Ty);  qmagRef = hypot(ref.Tx, ref.Ty);
errAll = max(abs([Tx(:) - ref.Tx(:); Ty(:) - ref.Ty(:)]));
[eLine, kLine] = max(abs(Tx(:) - ref.Tx(:)));  [iL, jL] = ind2sub([3 3], kLine);
fprintf('headline vs reference: max error %.2f K/m; largest dT/dx error at node (%d,%d): %.2f K/m (%.1f%% of the reference)\n', ...
        errAll, iL, jL, eLine, 100*eLine/abs(ref.Tx(kLine)));
fprintf('headline |q|/k vs reference: largest relative deviation %.1f%%\n', 100*max(abs(qmag(:)./qmagRef(:) - 1)));
fprintf('largest data error of the headline schemes: %.2f K/m\n', ...
        max(arrayfun(@(g) g.dataErr(strcmp(g.cand, g.best)), G)));
fprintf('headline q/k = -grad T, rows y = 3..1:\n  q_x/k\n');  disp(flipud(-Tx'));
fprintf('  q_y/k\n');  disp(flipud(-Ty'));
save(fullfile(resDir, 'b3_compare.mat'), 'G', 'Tx', 'Ty', 'schX', 'schY', 'qmag', 'qmagRef', 'errAll');

% ---- Figure: temperature along the row y = 2 and the slopes the schemes give at x = 1 and x = 2
f = figure('Visible', 'off');  hold on;
n = size(B2.Tfine, 1) - 1;  xf = linspace(0, P.L, n+1);  Jf = round(2/P.L*n) + 1;
plot(xf, B2.Tfine(:, Jf), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.6);
xg = 0:4;  Tg = [P.Tb, B1.Tcub(:,2)', P.Tb];
plot(xg, Tg, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
plot(P.xLine, P.Ts, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
st = {'C2', 'OS1', 'OS2', 'NU'};  ls = {':', '--', '-.', '-'};  col = lines(4);  hs = gobjects(1, 5);
for q = 1:numel(st)
    for i = 1:2
        s = B1.Dcub.(st{q}).Tx(i, 2);  if isnan(s), continue, end
        xs = i + [-0.3 0.3];  hs(q) = plot(xs, Tg(i+1) + s*(xs - i), ls{q}, 'Color', col(q,:), 'LineWidth', 1.3);
    end
end
for i = 1:2
    xs = i + [-0.3 0.3];  Ti = interp1(xf, B2.Tfine(:, Jf), i);
    hs(5) = plot(xs, Ti + ref.Tx(i,2)*(xs - i), 'k-', 'LineWidth', 2);
end
xline(P.xLine, 'r:');  box on;  grid on;  set(gca, 'FontSize', 8);  xlim([0 4]);  ylim([190 360]);
xlabel('$x$ [m]', 'Interpreter', 'latex');  ylabel('$T$ at $y=2$ [K]', 'Interpreter', 'latex');
legend(hs, {'N0 (central across line)', 'one-sided $O(h)$', 'one-sided $O(h^2)$', 'non-uniform NU', 'reference'}, ...
       'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 6.5);
saveFigB(f, 'b_schemes', 12, 7);

% ---- Figure: q/k arrows on the converged temperature map
f = figure('Visible', 'off');
imagesc(xf, xf, B2.Tfine', [200 350]);  set(gca, 'YDir', 'normal');  hold on;  colormap(turbo);
cb = colorbar;  cb.Label.String = 'T [K]';
plot([P.xLine P.xLine], [P.yLo P.yHi], 'k-', 'LineWidth', 2);
[Xn, Yn] = ndgrid(1:3, 1:3);  sc = 0.55 / max([qmag(:); qmagRef(:)]);
quiver(Xn(:), Yn(:), -sc*ref.Tx(:), -sc*ref.Ty(:), 0, 'w', 'LineWidth', 1.6, 'MaxHeadSize', 0.6);
quiver(Xn(:), Yn(:), -sc*Tx(:),     -sc*Ty(:),     0, 'k', 'LineWidth', 1.0, 'MaxHeadSize', 0.6);
plot(Xn(:), Yn(:), 'wo', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 3);
axis equal tight;  xlabel('x [m]');  ylabel('y [m]');  set(gca, 'FontSize', 8);
saveFigB(f, 'b_quiver', 9, 7);
