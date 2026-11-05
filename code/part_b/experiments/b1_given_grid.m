% B1_GIVEN_GRID  Gradient at the nine nodes from the given data (h = 1) by every scheme.
%   Reads part (a)'s nodal temperatures (does not re-run part a). Runs on T_cubic_h1, the answer
%   of part (a), and on T_reference, to measure how much part (a)'s error moves grad T.
here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, '..', 'results');  if ~exist(resDir, 'dir'), mkdir(resDir); end
A = readmatrix(fullfile(here, '..', '..', 'part_a', 'results', 'part_a_temperature.csv'));   % i,j,x,y,T_cubic_h1,T_sw_h1,T_reference
Tcub = zeros(3);  Tref9 = zeros(3);
for r = 1:size(A, 1)
    Tcub(A(r,1), A(r,2))  = A(r,5);
    Tref9(A(r,1), A(r,2)) = A(r,7);
end
Dcub = gradGiven(Tcub);
Dref = gradGiven(Tref9);

schemes = fieldnames(Dcub);
fprintf('\nGradient from the given grid (T_cubic_h1) [K/m], rows y = 3..1, cols x = 1..3\n');
for s = 1:numel(schemes)
    fprintf('%-4s dT/dx\n', schemes{s});  disp(flipud(Dcub.(schemes{s}).Tx'));
    fprintf('%-4s dT/dy\n', schemes{s});  disp(flipud(Dcub.(schemes{s}).Ty'));
end
fprintf('Richardson equals C4: max |R - C4| = %.2e\n', ...
    max(abs([Dcub.R.Tx(:) - Dcub.C4.Tx(:); Dcub.R.Ty(:) - Dcub.C4.Ty(:)]), [], 'omitnan'));
save(fullfile(resDir, 'b1_given_grid.mat'), 'Tcub', 'Tref9', 'Dcub', 'Dref');
