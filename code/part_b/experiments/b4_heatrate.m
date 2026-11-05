% B4_HEATRATE  Total heat rate per unit conductivity, Q/k, through the four edges (input for part c).
%   Given grid (h = 1): edge formula O(h^2) (2-point where the 3-point stencil would cross the line)
%   and O(h) everywhere, on T_cubic_h1 and on T_reference.  Refined fields: from b2_reference.
P = plateProblem();
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
B1 = load(fullfile(resDir, 'b1_given_grid.mat'));
B2 = load(fullfile(resDir, 'b2_reference.mat'));
grid5 = @(T9) [P.Tb*ones(1,5); P.Tb*ones(3,1), T9, P.Tb*ones(3,1); P.Tb*ones(1,5)];   % add the edges
H = struct();
[H.cub2, H.cub2E] = edgeFlux(grid5(B1.Tcub),  P.h, 2);
[H.ref2, H.ref2E] = edgeFlux(grid5(B1.Tref9), P.h, 2);
[H.cub1, H.cub1E] = edgeFlux(grid5(B1.Tcub),  P.h, 1);
[H.ref1, H.ref1E] = edgeFlux(grid5(B1.Tref9), P.h, 1);
fprintf('given grid, O(h^2) edge formula: Q/k = %.2f (T_cubic_h1), %.2f (T_reference)\n', H.cub2, H.ref2);
fprintf('given grid, O(h)   edge formula: Q/k = %.2f (T_cubic_h1), %.2f (T_reference)\n', H.cub1, H.ref1);
fprintf('per edge [east west north south]: O(h^2) %s   O(h) %s   reference %s\n', ...
        mat2str(H.cub2E, 4), mat2str(H.cub1E, 4), mat2str(B2.QEref, 4));
% headline edge formula: same rule as for the gradient (b3_compare) - the one closer to the reference
if abs(H.cub1 - B2.Qref) < abs(H.cub2 - B2.Qref), H.best = 1;  H.Qbest = H.cub1;  H.QbestE = H.cub1E;
else,                                             H.best = 2;  H.Qbest = H.cub2;  H.QbestE = H.cub2E; end
fprintf('refined fields: Q*/k = %.2f K m;  given grid, chosen O(h^%d) edge formula: %.2f (%+.1f%%)\n', ...
        B2.Qref, H.best, H.Qbest, 100*(H.Qbest/B2.Qref - 1));
save(fullfile(resDir, 'b4_heatrate.mat'), '-struct', 'H');

% ---- interface file for part (c):  h = (Q/k) * integral of k(t) dt
fid = fopen(fullfile(resDir, 'part_b_heatrate.csv'), 'w');
fprintf(fid, 'name,value,unit,description\n');
ord = {'O(h)', 'O(h^2)'};
fprintf(fid, 'Q_over_k_h1,%.4f,K m,given grid h=1 (T_cubic_h1) - %s edge formula (closer to the reference) + trapezoid\n', H.Qbest, ord{H.best});
fprintf(fid, 'Q_over_k_ref,%.4f,K m,refined fields N=64-128-256 + Richardson (best estimate)\n', B2.Qref);
fprintf(fid, 'Q_over_k_h1_O1,%.4f,K m,given grid - 2-point O(h) edge formula\n', H.cub1);
fprintf(fid, 'Q_over_k_h1_O2,%.4f,K m,given grid - 3-point O(h^2) edge formula (2-point where it would cross the line)\n', H.cub2);
names = {'east', 'west', 'north', 'south'};
for e = 1:4, fprintf(fid, 'Q_over_k_h1_%s,%.4f,K m,given grid (chosen formula) - %s edge only\n', names{e}, H.QbestE(e), names{e}); end
for e = 1:4, fprintf(fid, 'Q_over_k_ref_%s,%.4f,K m,refined fields - %s edge only\n', names{e}, B2.QEref(e), names{e}); end
fclose(fid);
