function saveFig(fig, name, w, h)
% SAVEFIG  Export figure FIG as vector PDF into report/figures/NAME.pdf, size w x h [cm].
dest = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'report', 'figures');
if ~exist(dest, 'dir'), mkdir(dest); end
theme(fig, 'light');                         % MATLAB's default theme may be dark; reports need white
set(fig, 'Units', 'centimeters', 'Position', [2 2 w h]);
exportgraphics(fig, fullfile(dest, [name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'white');
close(fig);
end
