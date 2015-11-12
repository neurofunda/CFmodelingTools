function view = ccmLoadDefault(view)
% ccmLoadDefault: load the default maps from a cortico-cortical model.
%
% 2009 KVH: adapted from rmLoadDefault.

if notDefined('view'), view = getCurView; end

view = ccmLoadModel(view, 1, 'varexplained', 'co');
view = ccmLoadModel(view, 1, 'sigmamajor', 'amp');

view = setCothresh(view, .15);
view = setClipMode(view, 'amp', ...
    [0 max(view.ccm.params.analysis.sigmas)]);

view.ui.coMode = setColormap(view.ui.coMode, 'hotCmap');
view.ui.ampMode = setColormap(view.ui.ampMode, 'jetCmap');

updateGlobal(view);
view  = refreshScreen(view);

return
