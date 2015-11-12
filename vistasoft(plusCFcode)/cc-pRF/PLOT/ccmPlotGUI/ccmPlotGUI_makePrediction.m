function [prediction, RF, rfParams, varexp,CFtSeries] = ccmPlotGUI_makePrediction(M, coords, voxel)
% Refresh the cortico-cortical model plotting GUI.
%
% 2011 KVH: wrote it

if notDefined('voxel')
    voxel = get(M.ui.voxel.sliderHandle, 'Value');
end
    
% Get model and info
model = M.model{M.modelNum};

% get variance explained
varexp = ccmCoordsGet(M.viewType, model, 'varexp', coords);

% varexp of zero indicates that there is no pRF: all the params are 0. If
% we go through the code below, it will error (NaNs get introduced). Don't
% bother -- we know the return values are empty:
if varexp == 0
	nTimePoints = size(M.params.analysis.allstimimages, 1);
	prediction = zeros(nTimePoints, 1);
    RF = zeros(numel(M.params.roi.coordsIndex), 1);
	rfParams = zeros(1,4);
	return
end

% get RF parameters from the model
rfParams = ccmPlotGUI_getRFParams(model, M.viewType, coords);

%% make predictions for each RF 

[trends, nt, dcid] = rmMakeTrends(M.params,0); 
trendBetas = pinv(trends)*M.params.analysis.allstimimages; 
allstimimages = M.params.analysis.allstimimages - trends*trendBetas; 


good_voxels = find(~isnan(allstimimages(1,:)));

for i = good_voxels     
    if rfParams(1) == M.sCoords(1,i) && ...
       rfParams(2) == M.sCoords(2,i) && ...
       rfParams(3) == M.sCoords(3,i) 
        X = M.distances(:,i);
        RF = single(exp(-1.*((X.^2)./(2.*rfParams(4).^2))));
        break;
    end
end

pred = allstimimages(:,good_voxels)*RF; 
%CFtSeries = allstimimages(:,M.modelNum);
CFtSeries = allstimimages;

% figure,plot(CFtSeries);
% figure,plot(pred);
% CFtSeries = allstimimages(M.modelNum);
% figure,plot(CFtSeries);

%% Compute final predicted time series (and get beta values)
beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel); 
prediction = [pred trends(:,dcid)] * beta; 


% recompute variance explained (because gridfit was with lscov)
rss = sum((M.tSeries(:,voxel)-prediction).^2);
rawrss = sum(M.tSeries(:,voxel).^2);

varexp = max(1 - rss./rawrss, 0); 

return
