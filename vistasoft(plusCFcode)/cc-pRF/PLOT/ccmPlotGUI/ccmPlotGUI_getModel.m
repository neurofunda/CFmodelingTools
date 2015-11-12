function M = ccmPlotGUI_getModel(view, roi)
% Retrieves the model data for the specified retinotopy model, in a struct
% M.
%
% 2011 KVH: adapted from rmPlotGUI_getModel.

% initialize the M struct with obvious fields

M.roiList=viewGet(view, 'roinames');
M.sourceROI = find(strcmp(roiList,view.ccm.params.sourceROI));
M.sCoords    = view.ROIs(M.sourceROI).coords;
M.roi        = roi;
M.model      = view.ccm.models;
M.modelNum   = 1;
M.params     = view.ccm.params;
M.dataType   = getDataTypeName(view);
M.viewType   = view.viewType;
M.distances  = ccmVoxToVoxDist(view.ROIs(M.sourceROI), view, view.mmPerVox);
M.prevVoxel  = 1;

% get time series and roi-coords
[M.tSeries, M.coords, M.params] = rmLoadTSeries(view, M.params, roi);


if M.params.analysis.fc
    M.tSeries = lowpass(M.tSeries,M.params);
end

% detrend
trends  = rmMakeTrends(M.params);
b = pinv(trends)*M.tSeries;
M.tSeries = M.tSeries - trends*b;

% make x-axis for time series plot
nFramesInd = [M.params.stim(:).nFrames]./[M.params.stim(:).nUniqueRep];
TR = [M.params.stim(:).framePeriod];
nScans = length(M.params.stim);

x = (0:nFramesInd(1)-1)'.*TR(1);
for n=2:nScans,
    timend = x(end);
    x = [x; ((0:nFramesInd(n)-1)'.*TR(n))+timend];  %#ok<AGROW>
end;
M.sepx = cumsum(nFramesInd-1) .* TR;
M.x    = x;

return;

function tSeries = lowpass(tSeries, params)
TR = params.stim.framePeriod;
fc = params.analysis.fc;
fs = 1/TR;
[b, a] = butter(10,fc/(fs/2),'low');
tSeries = filtfilt(b,a,tSeries);
return