% RUN_PART_B  Reproduce every number, table and figure of part (b).  Run from any folder.
%   Needs part (a)'s results/part_a_temperature.csv and its function solveFine (not re-run here).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'shared'));
addpath(fullfile(here, '..', 'part_a', 'core'), fullfile(here, '..', 'part_a', 'solvers'));
addpath(genpath(here));
resDir = fullfile(here, 'results');  if ~exist(resDir, 'dir'), mkdir(resDir); end
logFile = fullfile(resDir, 'run_log.txt');  if exist(logFile, 'file'), delete(logFile); end
diary(logFile);                           % console output is saved next to the results
tic;
b1_given_grid;
b2_reference;
b3_compare;
b4_heatrate;
test_part_b;                              % writes test_order.csv used by the tables
b5_tables;
fprintf('\nTotal run time %.1f s\n', toc);
diary off;
