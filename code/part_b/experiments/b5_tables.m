% B5_TABLES  Turn the saved part (b) results into LaTeX tabulars (tables/b_*.tex) and write the
%   nodal interface file results/part_b_flux.csv, so that no number in part (b) is typed by hand.
P = plateProblem();
resDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
tabDir = reportDir('tables');
B1 = load(fullfile(resDir, 'b1_given_grid.mat'));
B2 = load(fullfile(resDir, 'b2_reference.mat'));
B3 = load(fullfile(resDir, 'b3_compare.mat'));
B4 = load(fullfile(resDir, 'b4_heatrate.mat'));
label = struct('C2', 'C2', 'C4', 'C4', 'R', 'R', 'OS1', 'OS1', 'OS2', 'OS2', 'NU', 'NU');

% ---- Answer table: q/k = -grad T laid out like the grid, headline in bold, reference in brackets
fid = openTab(tabDir, 'b_answer', '@{}l c|ccc@{}', ' & $y \backslash x$ & $x=1$ & $x=2$ & $x=3$');
blocks = {'$q_x/k$', -B3.Tx, -B2.GXref; '$q_y/k$', -B3.Ty, -B2.GYref; '$|\vec q|/k$', B3.qmag, B3.qmagRef};
for b = 1:3
    if b > 1, fprintf(fid, '\\midrule\n'); end
    for j = 3:-1:1
        lab = '';  if j == 3, lab = blocks{b,1}; end
        fprintf(fid, '%s & $y=%d$', lab, j);
        for i = 1:3, fprintf(fid, ' & $\\mathbf{%.1f}$ [$%.1f$]', blocks{b,2}(i,j), blocks{b,3}(i,j)); end
        fprintf(fid, ' \\\\\n');
    end
end
closeTab(fid);

% ---- Scheme table: dT/dx at the four nodes beside the line, every scheme, with % deviation
st = {'C2', 'OS1', 'OS2', 'NU'};
fid = openTab(tabDir, 'b_schemes', '@{}l rrrr r@{}', ...
    'Node & N0 (central) & OS1, $O(h)$ & OS2, $O(h^2)$ & NU, $O(h^2)$ & Reference');
for nd = [1 2; 1 3; 2 2; 2 3]'
    i = nd(1);  j = nd(2);  r = B2.GXref(i,j);
    fprintf(fid, '$(%d,%d)$', i, j);
    for q = 1:numel(st)
        v = B1.Dcub.(st{q}).Tx(i,j);
        if isnan(v), fprintf(fid, ' & n/a');
        else, fprintf(fid, ' & $%.1f$ ($%+.0f\\%%$)', v, 100*(v/r - 1)); end
    end
    fprintf(fid, ' & $%.1f$ \\\\\n', r);
end
closeTab(fid);

% ---- Refinement table: gradient at representative nodes and Q/k for N = 32..256, order, Richardson
fid = openTab(tabDir, 'b_reference', '@{}r r rrrr r@{}', ...
    '$N$ & $h$ [m] & $T_x(1,2)$ & $T_x(2,2)$ & $T_y(1,3)$ & $T_y(2,2)$ & $Q/k$ [K\,m]');
for k = 1:numel(B2.Ns)
    fprintf(fid, '$%d$ & %.4f & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\\n', B2.Ns(k), P.L/B2.Ns(k), ...
        B2.GX(k,1,2), B2.GX(k,2,2), B2.GY(k,1,3), B2.GY(k,2,2), B2.Qk(k));
end
fprintf(fid, '\\midrule\nobserved $p$ & & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
    B2.pX(1,2), B2.pX(2,2), B2.pY(1,3), B2.pY(2,2), B2.pQ);
fprintf(fid, 'Richardson & & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\\n', ...
    B2.GXref(1,2), B2.GXref(2,2), B2.GYref(1,3), B2.GYref(2,2), B2.Qref);
closeTab(fid);

% ---- Error table: per group, every candidate; chosen scheme in bold
gname = {'$T_x$ at $(1,2),(1,3)$', '$T_x$ at $(2,2),(2,3)$', '$T_x$ at $(2,1)$', ...
         '$T_x$ at $(1,1),(3,j)$', '$T_y$ on row $y=2$', '$T_y$ on rows $y=1,3$'};
order = struct('C2', '2', 'C4', '4', 'OS1', '1', 'OS2', '2', 'NU', '2');
fid = openTab(tabDir, 'b_errors', '@{}l l c rrr@{}', ...
    'Nodes & Scheme & Order & Total & Data & Differentiation');
for g = 1:numel(B3.G)
    if g > 1, fprintf(fid, '\\midrule\n'); end
    for c = 1:numel(B3.G(g).cand)
        s = B3.G(g).cand{c};  lab = '';  if c == 1, lab = gname{g}; end
        nm = label.(s);  if strcmp(s, 'C2') && g <= 2, nm = 'N0 (C2 across line)'; end
        fmt = '%s & %s & %s & %.1f & %.1f & %.1f \\\\\n';
        if strcmp(s, B3.G(g).best), fmt = '%s & \\textbf{%s} & %s & \\textbf{%.1f} & %.1f & %.1f \\\\\n'; end
        fprintf(fid, fmt, lab, nm, order.(s), B3.G(g).err(c), B3.G(g).dataErr(c), B3.G(g).diffErr(c));
    end
end
closeTab(fid);

% ---- Heat-rate table
fid = openTab(tabDir, 'b_heatrate', '@{}l rrrrr r@{}', ...
    'Data and edge formula & $Q/k$ & East & West & North & South & vs.\ $Q^*$');
rows = {'$h=1$, \texttt{T\_cubic\_h1}, $O(h^2)$', B4.cub2, B4.cub2E; ...
        '$h=1$, \texttt{T\_cubic\_h1}, $O(h)$',   B4.cub1, B4.cub1E; ...
        '$h=1$, \texttt{T\_reference}, $O(h^2)$', B4.ref2, B4.ref2E; ...
        '$h=1$, \texttt{T\_reference}, $O(h)$',   B4.ref1, B4.ref1E};
rows{3 - B4.best, 1} = ['\textbf{' rows{3 - B4.best, 1} '}'];   % row of the chosen formula in bold
for k = 1:numel(B2.Ns)
    rows(end+1,:) = {sprintf('refined $N=%d$, $O(h^2)$', B2.Ns(k)), B2.Qk(k), B2.QE(k,:)}; %#ok<SAGROW>
end
rows(end+1,:) = {'Richardson $Q^*/k$', B2.Qref, B2.QEref};
for k = 1:size(rows, 1)
    if k == 5 || k == size(rows, 1), fprintf(fid, '\\midrule\n'); end
    fprintf(fid, '%s & %.1f & %.1f & %.1f & %.1f & %.1f & $%+.1f\\%%$ \\\\\n', rows{k,1}, rows{k,2}, rows{k,3}, ...
            100*(rows{k,2}/B2.Qref - 1));
end
closeTab(fid);

% ---- Order-of-accuracy test table (from test_part_b)
C = readcell(fullfile(resDir, 'test_order.csv'));
fid = openTab(tabDir, 'b_test', '@{}l c r r r@{}', 'Scheme & Expected order & Error at $h=0.4$ & Error at $h=0.025$ & Observed order');
for k = 2:size(C, 1)
    fprintf(fid, '%s & %d & $%s$ & $%s$ & %.2f \\\\\n', C{k,1}, C{k,2}, sci(C{k,3}), sci(C{k,4}), C{k,5});
end
closeTab(fid);

% ---- nodal interface file
fid = fopen(fullfile(resDir, 'part_b_flux.csv'), 'w');
fprintf(fid, 'i,j,x,y,Tx,Ty,scheme_x,scheme_y,qx_over_k,qy_over_k,qmag_over_k,Tx_ref,Ty_ref,qmag_ref_over_k\n');
for j = 1:3, for i = 1:3
    fprintf(fid, '%d,%d,%d,%d,%.4f,%.4f,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', i, j, i, j, B3.Tx(i,j), B3.Ty(i,j), ...
        B3.schX{i,j}, B3.schY{i,j}, -B3.Tx(i,j), -B3.Ty(i,j), B3.qmag(i,j), B2.GXref(i,j), B2.GYref(i,j), B3.qmagRef(i,j));
end, end
fclose(fid);

function fid = openTab(tabDir, name, spec, header)
% Open NAME.tex and write the tabular preamble (column spec, top rule, header row, mid rule).
fid = fopen(fullfile(tabDir, [name '.tex']), 'w');
fprintf(fid, '\\begin{tabular}{%s}\n\\toprule\n%s \\\\\n\\midrule\n', spec, header);
end

function closeTab(fid)
% Write the bottom rule and close the tabular.
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');  fclose(fid);
end

function s = sci(v)
% LaTeX scientific notation with two decimals, e.g. 2.98\times10^{-8}.
e = floor(log10(abs(v)));  s = sprintf('%.2f\\times10^{%d}', v/10^e, e);
end
