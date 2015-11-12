function [data, params, coords] = ccmLoadData_mod(view, params, roi)
% ccmLoadData - load and preprocess time series
%
% 2009: KVH adapted from rmLoadData.

data = [];

% place datasets behind each other.
for ds = 1:numel(params.stim),
    [tSeries coords params] = ccmLoadDataROI(view, params, ds, roi);
    
     if params.win == 1
        
        X = [1:1:length(tSeries(:,1))]';
        win = taper(X,'TaperPar',0.5);
%         figure,plot(win);
%         figure,subplot(211);plot(tSeries);
        for i = 1:length(tSeries(1,:));
            tSeries(:,i) = tSeries(:,i).*win;
        end
%         subplot(212);plot(tSeries);
     end
    
     
    
    if isempty(data),
        dii.end = cumsum([params.stim(:).nFrames]./[params.stim(:).nUniqueRep]);
        dii.start = [1 dii.end(1:end-1)+1];
        data = zeros(dii.end(end), size(tSeries ,2));
        
    end;
    
   
    data(dii.start(ds):dii.end(ds),:) = tSeries;

    
end;

return;


function data=raw2pc(data)
% convert to percent signal change
dc = ones(size(data,1),1)*mean(data);
data = ((data./dc) - 1) .*100;
return;


function [tSeries coords params] = ccmLoadDataROI(view, params, ds, r)

% get ROI coords
coords = view.ROIs(r).coords;

% index into view's data
[coordsIndex coords] = roiIndices(view, coords);

% store roi info
if r == view.selectedROI,
    params = ccmSet(params,'roiName',view.ROIs(r).name);
    params = ccmSet(params,'roiCoords',coords);
    params = ccmSet(params,'roiIndex',coordsIndex);
end


tSeries  = loadtSeries(view, ds);

roiIndex = coordsIndex;

if params.surr==1
    
    %% iAAFT surrogates
    
    tSeries = tSeries(:,roiIndex(:));
    for i=1:length(tSeries)
        tSeries(:,i) = surrogate(tSeries(:,i),1); % iAAFF surrogates
    end
    
    coords = roiIndex(:);
    
else
    %% normal step  tSeries = tSeries(:,roiIndex(:)); coords = roiIndex(:);
   
    tSeries = tSeries(:,roiIndex(:));
    coords = roiIndex(:);
    
    
end

% only convert to percent change if the flag is set
if params.analysis.calcPC
    tSeries = raw2pc(tSeries);
    %% Fig   unfiltered (1)
    %figure, plot(tSeries);
end

% apply low-pass filtering
if params.analysis.fc
    tSeries = ccmLowPass(tSeries,params);
    %% Fig  LP filtered (2)
    %figure, plot(tSeries);
    
end

% remove global signals
if params.analysis.gsr
    tSeries = ccmGlobalSignalRegress(view,ds,tSeries,1);
    %% Fig  GS regressed (not applied)
    %figure, plot(tSeries);
end

return

