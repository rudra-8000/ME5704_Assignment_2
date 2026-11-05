function d = reportDir(kind)
% REPORTDIR  Folder of the LaTeX report for KIND = 'figures' or 'tables' (created if missing).
%   Overleaf project: code/ sits next to sections/, figures/ and tables/.
%   Repository:       code/ and report/ are siblings.
root = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
if exist(fullfile(root, 'sections'), 'dir'), d = fullfile(root, kind);
else, d = fullfile(root, 'report', kind); end
if ~exist(d, 'dir'), mkdir(d); end
end
