function params = ccmFit(view,params)
% ccmFit - brute force fitting of predictions based upon premade cc-pRFs
%
% 2009 KVH: wrote it.
% argument checks
if notDefined('view'),   error('Need view struct'); end;
if notDefined('params'), error('Need params'); end;

% For speed we do our computations in single precision, but we output in
% double (for compatibility).
params.analysis.allstimimages = single(params.analysis.allstimimages);

% make trends to fit with the model (discrete cosine set)
[trends, ntrends, dcid] = rmMakeTrends(params);
trends = rmDecimate(trends,params.analysis.decimate);
trends = single(trends);

% make all predictions first
allstimimages = rmDecimate(params.analysis.allstimimages, params.analysis.decimate);

% detrend prediction
trendBetas = pinv(trends)*allstimimages;
allstimimages = allstimimages - trends*trendBetas;

% compute the distance between each pair of voxels in the source ROI
Distances = ccmVoxToVoxDist(view.ROIs(params.sourceROI), view, view.mmPerVox);

% remove voxels without signal (this shouldn't happen, but can happen when
% importing scans for other sessions).
good_voxels = find(~isnan(params.analysis.allstimimages(1,:)));

% % if insufficient memory we could try decreasing the number of cc-pRF sizes
% mem = memory;
% if numel(good_voxels)*numel(good_voxels)*length(params.analysis.sigmas)*4 > mem.MaxPossibleArrayBytes
%     params.analysis.sigmas = linspace(min(params.analysis.sigmas),max(params.analysis.sigmas),numel(params.analysis.sigmas)/2);
% end
%params.analysis.sigmas = linspace(1,20,numel(params.analysis.sigmas));

% for speed we predefine some variables
S = single(repmat(params.analysis.sigmas(:)',numel(good_voxels),1));
RF = single(zeros(numel(good_voxels),numel(good_voxels)*length(params.analysis.sigmas)));
x = single(zeros(1,numel(good_voxels)*length(params.analysis.sigmas)));
y = x; z = x; s = x;

% user feedback
fprintf(1,'[%s]:Computing %s model predictions [%s]: ... ', mfilename, ...
    num2str(numel(params.analysis.centers)*numel(params.analysis.sigmas)), ...
    view.ROIs(params.sourceROI).name); drawnow;tic;

% define the cc-PRF models
for i = good_voxels
    ii = (1:length(params.analysis.sigmas))+length(params.analysis.sigmas)*(i-1);
    x(ii) = single(repmat(view.ROIs(params.sourceROI).coords(1,i),1,numel(params.analysis.sigmas)));
    y(ii) = single(repmat(view.ROIs(params.sourceROI).coords(2,i),1,numel(params.analysis.sigmas)));
    z(ii) = single(repmat(view.ROIs(params.sourceROI).coords(3,i),1,numel(params.analysis.sigmas)));
    s(ii) = single(params.analysis.sigmas);
    X = single(repmat(Distances(good_voxels,i),1,numel(params.analysis.sigmas)));
    RF(:,ii) = single(exp(-1.*((X.^2)./(2.*S.^2))));
end

% user feedback
fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now)); drawnow;

% generate the cc-pRF model predictions
prediction = allstimimages(:,good_voxels)*RF;

% free some very large memory chunks
clear Distances RF

% add the cc-pRF model parameters to the output struct
params.analysis.x = single(x);
params.analysis.y = single(y);
params.analysis.z = single(z);
params.analysis.sigmaMajor = single(s);
% free some more memory
clear x y z s

% load the data
[data, params] = ccmLoadData(view, params, params.targetROI);


% for speed convert to single and remove NaNs (should not be there anyway!)
data(isnan(data)) = 0;
data = rmDecimate(data,params.analysis.decimate);
data = single(data);



% remove trends from data so they do not count in the percent variance
% explained calculation later.
trendBetas = pinv(trends)*data;
data = data - trends*trendBetas;

% compute rss raw data for variance computation later
rssdata = sum(data.^2);

% initiate model
fprintf(1,'[%s]:Number of voxels: %d.\n',mfilename,size(data,2));drawnow;

%% initialize model
model = initiateModel(params, size(data,2), ntrends);

for mm = 1:numel(model),
    model{mm} = ccmSet(model{mm},'roiCoords',ccmGet(params,'roiCoords'));
    model{mm} = ccmSet(model{mm},'roiIndex',ccmGet(params,'roiIndex'));
    model{mm} = ccmSet(model{mm},'roiName',ccmGet(params,'roiName'));
end;

% put in number of data points.
for mm = 1:numel(model),
    model{mm} = ccmSet(model{mm},'npoints',size(data,1));
end;

% now we extract only the data from that slice and put it in a
% temporary structure that will be modified throughout.
s = ccmSliceGet(model);

% initiateModel fills the rss-field with Infs. We reset them here
% to a more data-driven maximum value of sum(data.^2)
for n=1:numel(s),
    s{n}.rawrss = rssdata;
end;

% fit different cc-pRF models

t.trends = trends(:,dcid);
t.dcid = dcid;
for n=1:numel(params.analysis.pRFmodel)
    s{n}=ccmFit_core(s{n},prediction,data,params,t);
end

% now put back the trends to the fits
for mm=1:numel(s),
    nB = size(s{mm}.b,1);
    s{mm}.b(nB-ntrends+1:end,:) = s{mm}.b(nB-ntrends+1:end,:)+trendBetas;
end

% now we put back the temporary data from that slice
model = ccmSliceSet(model,s);

% save
ccmSave(view,model,params);

return;


function model = initiateModel(params,d2,nt)
% make the model struct with ccmSet
fillwithzeros = zeros(1,d2);
fillwithinfs = ones(1,d2).*Inf;

% add a small number to sigmas because a cc-pRF with 0 sigma does not exist
smallnumber =  0.001 ;

% initiate all models
model = cell(numel(params.analysis.pRFmodel),1);
for n=1:numel(params.analysis.pRFmodel),
    model{n} = ccmSet;
    model{n} = ccmSet(model{n},'x0',fillwithzeros);
    model{n} = ccmSet(model{n},'y0',fillwithzeros);
    model{n} = ccmSet(model{n},'z0',fillwithzeros);
    model{n} = ccmSet(model{n},'s_major',fillwithzeros+smallnumber);
    model{n} = ccmSet(model{n},'rawrss',fillwithzeros);
    model{n} = ccmSet(model{n},'rss',fillwithinfs);
    model{n} = ccmSet(model{n},'ntrends',nt);
    model{n} = ccmSet(model{n},'b',zeros(1,d2,nt+1));
    model{n} = ccmSet(model{n},'desc','2D pRF fit (x,y,sigma, positive only)');
end;

return;
