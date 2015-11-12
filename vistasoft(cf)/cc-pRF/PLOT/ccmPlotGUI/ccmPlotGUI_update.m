function ccmPlotGUI_update(vw, M)
% Refresh the cortico-cortical model plotting GUI.
%
% 2011 KVH: adapted from ccmPlotGUI_update

if ~exist('M','var') || isempty(M), M = get(gcf, 'UserData'); end

% get current voxel from the GUI
voxel = get(M.ui.voxel.sliderHandle, 'Value');

% get needed params
coords = M.coords(:,voxel);
if isequal(M.roi.viewType, 'Gray')  % convert coords into an index
    coords = M.coords(voxel);
end

M.modelNum = get(M.ui.model, 'Value');

% compute RF, get tSeries for this voxel
[pred RF rfParams variance_explained] = ccmPlotGUI_makePrediction(M, coords, voxel);

% check if the prediction for this voxel is empty. If so, don't plot it:
if all( isnan(pred) | isinf(pred) )
	axes(M.ui.tsAxes);  cla;   %#ok<MAXES>
	AX = axis;
	text(AX(2) + .5*(AX(2)-AX(1)), AX(3) + .5*(AX(4)-AX(3)), ...
		 '(No data available for this voxel)', ...
		 'FontSize', 14, ... 
		 'HorizontalAlignment', 'center');
	return
end

% store the current prediction
M.prediction = pred;
M.RF = RF;

% (1) plot time series

axes(M.ui.tsAxes); cla; hold on; %#ok<MAXES>

hTSeries = plot(M.x, M.tSeries(:,voxel), 'k:', 'LineWidth', 1.5);
hFit = plot(M.x, pred(:, 1), 'b', 'LineWidth', 1.5);
hResidual = plot(M.x, M.tSeries(:,voxel)-pred(:,1), 'r:', 'LineWidth', 1);

allPlotted = [M.tSeries(:,voxel); pred(:,1); M.tSeries(:,voxel)-pred(:,1)];
axis([min(M.x) max(M.x) min(allPlotted) max(allPlotted)]);
h = axis;
for n=1:numel(M.sepx),
    plot([1 1].*M.sepx(n), [h(3) h(4)], 'k:', 'LineWidth', 2);
end;
xlabel('Time (sec)');
ylabel('BOLD signal change (%)');

% set the user data of each checkbox to point to the PLOT curves
% we just created. The callbacks for these checkboxes will then toggle
% the visibility of each curve.
set(M.ui.tsCheck, 'UserData', hTSeries);
set(M.ui.predCheck, 'UserData', hFit);
set(M.ui.resCheck, 'UserData', hResidual);
if get(M.ui.tsCheck, 'Value')==0, set(hTSeries, 'Visible', 'off'); end
if get(M.ui.predCheck, 'Value')==0, set(hFit, 'Visible', 'off'); end
if get(M.ui.resCheck, 'Value')==0, set(hResidual, 'Visible', 'off'); end

% (2) set text fields  

% print out percent variance (R^2) explained
txt = sprintf('Variance explained: %.2f%%', variance_explained*100);
set(M.ui.r2Text, 'String', txt);

% print out coords of voxel
txt = ['Coords: [' num2str(M.params.roi.coords(:,voxel)') ']'];
set(M.ui.coordsText, 'String', txt);

txt = sprintf('x=%.1f, y=%.1f, z=%.1f',rfParams(1), rfParams(2), rfParams(3));
set(M.ui.xyzText, 'String', txt);

txt = sprintf('sigma=%.1f%s', rfParams(4), 'mm');
set(M.ui.sigmaText, 'String', txt);

% from rmPlot: cache the data for the current time series plot as well
M.currTsData.x    = M.x;
M.currTsData.pred = pred;
M.currTsData.raw  = M.tSeries(:,voxel);
M.currTsData.sepx = M.sepx;
M.currTsData.pRF  = RF;

% store the updated M struct in the figure's user data
M.prevVoxel = voxel;  % also update the most recent voxel
set(M.fig, 'UserData', M);

% display the connective field on mesh if requested
if get(M.ui.showCheck, 'Value') == 1 && any(rfParams)
    ccmShowConnectiveField(vw, rfParams, M.distances);
end

return;
