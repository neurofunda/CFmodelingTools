function M = ccmPlotGUIsw_openFig(M,voxel)
% Open the ccmPlot GUI figure.
%
% 2011 KVH: adapted from rmPlotGUI_openfig.
% 2014 NG: modified it

% switch for dealing with java figures
%javaFigs = mrvJavaFeature;

bckg = [0.94 0.94 0.94];

% open the figure
figName = sprintf('ccmPlotGUI + time series analysis %s (%s)', M.roi.name, M.dataType);
M.fig = figure('Color', 'w', 'Name', figName, ...
    'Units', 'pixels', 'Position', [600 400 925 489]);

% create a set of axes for showing the time series / prediction / RSS
 %M.ui.tsAxes = axes('position',[.06  .57  .7  .35]);
% 
 %M.ui.cfAxes = axes('position',[.06  .08  .7  .35]);

M.ui.tsAxes = axes('position',[.06  .72  .7  .25]);

M.ui.corrAxes = axes('position',[.06  .4  .7  .2]);

M.ui.cfAxes = axes('position',[.06  .08  .7  .25]);

% M.ui.polarMapAxes = axes('position',[.6  .08  .2  .2]);

% create sliders / text controls
callback = 'ccmPlotGUIsw([],[],[],1);';

M.ui.model = uicontrol('Style', 'popupmenu', ...
    'Units', 'norm', 'Position', [.8 .87 .2 .1], ...
    'BackgroundColor',bckg, 'FontSize', 10, ...
    'Callback', callback, ...
    'String', 'CC-pRF model', 'Value', 1);


M.ui.voxel = mrvSlider([.8 .82 .2 .1], 'Voxel', ...
    'Range', [1 size(M.tSeries, 2)], 'IntFlag', 1, 'Value', voxel, ...
    'MaxLabelFlag', 1, 'FontSize', 10, ...
    'Color', bckg,'BackgroundColor', bckg, 'Callback', callback);

callback = ['opt = get(gcbo, ''Value'') + 1;  tmp = {''off'' ''on''}; ' ...
    'set( get(gcbo, ''UserData''), ''Visible'', tmp{opt} ); ' ...
    'clear opt tmp; '];


M.ui.tsCheck = uicontrol('Style', 'checkbox', ...
    'Units', 'norm', 'Position', [.8 .75 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'Callback', callback, 'Min', 0, 'Max', 1, ...
    'String', 'Time Series', 'Value', 1);

M.ui.predCheck = uicontrol('Style', 'checkbox', ...
    'Units', 'norm', 'Position', [.9 .75 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'Callback', callback, 'Min', 0, 'Max', 1, ...
    'String', 'Prediction', 'Value', 1);

M.ui.resCheck = uicontrol('Style', 'checkbox', ...
    'Units', 'norm', 'Position', [.8 .70 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'Callback', callback, 'Min', 0, 'Max', 1, ...
    'String', 'Residuals', 'Value', 0);

M.ui.showCFsignal = uicontrol('Style', 'checkbox', ...
    'Units', 'norm', 'Position', [.9 .70 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'Callback', callback, 'Min', 0, 'Max', 1, ...
    'String', 'CF Signal', 'Value', 1);

M.ui.showCheck = uicontrol('Style', 'checkbox', ...
    'Units', 'norm', 'Position', [.9 .65 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'Callback', callback, 'Min', 0, 'Max', 1, ...
    'String', 'Mesh', 'Value', 1);

M.ui.showMov = uicontrol('Style', 'checkbox', ...
    'Units', 'norm', 'Position', [.8 .65 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'Callback', callback, 'Min', 0, 'Max', 1, ...
    'String', 'Movie', 'Value', 0);


M.ui.mapCoverage = uicontrol('Style', 'checkbox', ...
    'Units', 'norm', 'Position', [.8 .6 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'Callback', callback, 'Min', 0, 'Max', 1, ...
    'String', 'Coverage', 'Value', 0);

M.ui.showConnections = uicontrol('Style', 'checkbox', ...
    'Units', 'norm', 'Position', [.9 .6 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'Callback', callback, 'Min', 0, 'Max', 1, ...
    'String', 'Connections', 'Value', 1);
% add text uicontrols which will contain outputs on the calculations

M.ui.r2Text = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .45 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');

M.ui.r2TextRef = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .425 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');


M.ui.coordsText = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .375 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');

M.ui.xyzText = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .35 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');

M.ui.eccText = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .3 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');

M.ui.eccTextRef = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .275 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');

M.ui.polText = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .225 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');

M.ui.polTextRef = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .2 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');

M.ui.sigmaText = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .15 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', ''); % 'Interpreter', 'tex',

M.ui.sigmaTextRef = uicontrol('Style', 'text', ...
    'HorizontalAlignment', 'left', ...
    'Units', 'norm', 'Position', [.8 .125 .2 .1], ...
    'BackgroundColor', bckg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'String', '');

% lastly, initialize the GUI graphics by performing an update:
set(M.fig, 'UserData', M);
%mrvJavaFeature(javaFigs);

return