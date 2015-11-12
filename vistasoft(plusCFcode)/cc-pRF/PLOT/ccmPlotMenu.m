function ccmMenu = ccmPlotMenu(h, view)
% Submenu for cortico-cortical pRF modeling
%
% 2012 KVH: wrote it.

ccmMenu = uimenu(h, 'Label', 'Cortico-Cortical Model', 'Separator', 'off');

% Plot the pRF linkage between two ROIs 
callback = sprintf('ccmPlotLinkage(%s);', view.name);
uimenu(ccmMenu, 'Label', 'Plot receptive field linkage', ...
    'Separator', 'off', 'CallBack', callback);

% Plot the cc-pRF size versus eccentricity for the current ROI
callback = sprintf('ccmPlotEccSigma(%s);', view.name);
uimenu(ccmMenu, 'Label', 'Plot cc-pRF size vs. pRF eccentricity (current ROI)', ...
    'Separator', 'on', 'CallBack', callback);

% Plot the cc-pRF size versus eccentricity for multiple ROIs
callback = sprintf('ccmPlotMultiEccSigma(%s);', view.name);
uimenu(ccmMenu, 'Label', 'Plot cc-pRF size vs. pRF eccentricity (selected ROIs)', ...
	'Separator', 'off', 'CallBack', callback);

% Plot the cc-pRF size versus laterality for the current ROI
callback = sprintf('ccmPlotLatSigma(%s);', view.name);
uimenu(ccmMenu, 'Label', 'Plot cc-pRF size vs. pRF laterality (current ROI)', ...
    'Separator', 'on', 'CallBack', callback);

% Plot the cc-pRF size versus laterality for multiple ROIs
callback = sprintf('ccmPlotMultiLatSigma(%s);', view.name);
uimenu(ccmMenu, 'Label', 'Plot cc-pRF size vs. pRF laterality (selected ROIs)', ...
	'Separator', 'off', 'CallBack', callback);

% Plot the sampling radius versus eccentricity for the current ROI
callback = sprintf('ccmPlotEccLatSigma(%s);', view.name);
uimenu(ccmMenu, 'Label', 'Plot the sampling radius vs. eccentricity (current ROI)', ...
    'Separator', 'on', 'CallBack', callback);

% Plot the sampling radius versus laterality for multiple ROIs
callback = sprintf('ccmPlotMultiEccLatSigma(%s);', view.name);
uimenu(ccmMenu, 'Label', 'Plot the sampling radius vs. eccentricity (selected ROIs)', ...
	'Separator', 'off', 'CallBack', callback);

% Visualize the model fit (all time points)
callback = ['ccmPlotGUI(', view.name, ',[]);'];
uimenu(ccmMenu, 'Label', 'Visualize the cortico-cortical model fit', ...
    'Separator', 'on', 'CallBack', callback);

return