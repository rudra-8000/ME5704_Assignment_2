% RUN_PART_A  Reproduce every number and figure of part (a).  Run from any folder.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'shared'));
addpath(genpath(here));
logFile = fullfile(here, 'results', 'run_log.txt');  if exist(logFile, 'file'), delete(logFile); end
diary(logFile);                           % console output is saved next to the results
rng(1);                                   % fixed seed -> reproducible Monte Carlo
tic;
a1_treatments;
a2_refinement;
a3_solvers;
a4_transient;
a5_montecarlo;
a6_tables;
test_part_a;
fprintf('\nTotal run time %.1f s\n', toc);
diary off;
