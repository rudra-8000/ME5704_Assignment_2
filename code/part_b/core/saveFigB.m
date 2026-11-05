function saveFigB(fig, name, w, h)
% SAVEFIGB  Same export as part (a)'s saveFig (vector PDF, w x h cm, white background), but into
%   the figures folder found by reportDir, so it works in the Overleaf layout as well.
theme(fig, 'light');
set(fig, 'Units', 'centimeters', 'Position', [2 2 w h]);
exportgraphics(fig, fullfile(reportDir('figures'), [name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'white');
close(fig);
end
