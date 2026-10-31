% A1_TREATMENTS  Six ways of using the given h = 1 data with the source line between nodes.
%   Solves each system with Gauss elimination, prints/saves the 3x3 fields, and draws
%   (i) the problem schematic and (ii) the six temperature maps.
P = plateProblem();
methods = {'sw','ghost1','ghost3','snapL','snapR','snapAll'};
titles  = {'Shortley-Weller','Ghost, linear','Ghost, cubic','Snap left pair','Snap right pair','Snap all four'};
resDir  = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
R = struct();
for m = 1:numel(methods)
    [A, b] = buildSystem(methods{m});
    R.(methods{m}) = reshape(gaussElim(A, b), 3, 3);     % R.method(i,j)
end
save(fullfile(resDir, 'a1_treatments.mat'), '-struct', 'R');

% ---- print + csv (long format: method,i,j,T)
fid = fopen(fullfile(resDir, 'a1_treatments.csv'), 'w');  fprintf(fid, 'method,i,j,T\n');
for m = 1:numel(methods)
    fprintf('\n%s (rows y=3..1, cols x=1..3)\n', titles{m});  disp(flipud(R.(methods{m})'));
    for j = 1:3, for i = 1:3
        fprintf(fid, '%s,%d,%d,%.4f\n', methods{m}, i, j, R.(methods{m})(i,j));
    end, end
end
fclose(fid);

% ---- Figure: schematic of the plate, grid and source line
f = figure('Visible','off');  hold on;  axis equal;  axis([-0.3 4.3 -0.3 4.6]);
rectangle('Position', [0 0 4 4], 'LineWidth', 1.5);
for i = 1:3, for j = 1:3
    plot(i, j, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    text(i-0.15, j+0.22, sprintf('(%d,%d)', i, j), 'FontSize', 8);
end, end
plot([P.xLine P.xLine], [P.yLo P.yHi], 'r-', 'LineWidth', 2.5);
plot([1 2], [2.5 2.5], 'b-');
text(1.12, 2.38, '0.5', 'Color', 'b', 'FontSize', 8);  text(1.67, 2.38, '0.5', 'Color', 'b', 'FontSize', 8);
text(1.55, 1.62, '$T_s = 350$ K', 'Color', 'r', 'Interpreter', 'latex', 'FontSize', 9);
text(0.6, 4.25, '$T = 200$ K on all edges', 'Interpreter', 'latex', 'FontSize', 9);
xlabel('x'); ylabel('y'); box on;  set(gca, 'FontSize', 9);
saveFig(f, 'a_schematic', 8, 7.5);

% ---- Figure: six 3x3 temperature maps
f = figure('Visible','off');
for m = 1:numel(methods)
    subplot(2, 3, m);  imagesc(1:3, 1:3, R.(methods{m})', [200 350]);  set(gca, 'YDir', 'normal');
    for i = 1:3, for j = 1:3
        text(i, j, sprintf('%.0f', R.(methods{m})(i,j)), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', txtColor(R.(methods{m})(i,j)));
    end, end
    title(titles{m}, 'FontSize', 8);  axis equal tight;  set(gca, 'XTick', 1:3, 'YTick', 1:3, 'FontSize', 8);
end
colormap(turbo);  cb = colorbar('Position', [0.93 0.2 0.015 0.6]);  cb.Label.String = 'T [K]';
saveFig(f, 'a_treatments', 15, 8);

function c = txtColor(T)
% White text on the dark ends of the colour map, black in the middle.
c = 'k';  if T < 240 || T > 335, c = 'w'; end
end
