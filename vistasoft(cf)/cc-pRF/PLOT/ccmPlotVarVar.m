function ccmPlotVarVar(vw)
% Compares cc-pRF varexp with pRF varexp
%
% KVH 2012: wrote it 

% cortico-cortical model must be loaded
try
	ccm = vw.ccm.models{vw.ccm.modelNum};
catch %#ok<CTCH>
    ccm = [];
end;

if isempty(ccm)
    fprintf(1,'[%s]: please load a cortico-cortical model first\n', ...
        mfilename); drawnow; return;
end

% select retModel
rmFile = getPathStrDialog(dataDir(vw),...
    'Choose pRF model file name', 'retModel*.mat');
tmp = load(rmFile);
rm = tmp.model{1};

% check ROIs
roiList=viewGet(vw, 'roinames');
if isempty(roiList)
    fprintf(1,'[%s]: please load and select at least one ROI for plotting\n', ...
        mfilename); drawnow; return;
end
selectedROIs = find(buttondlg('ROIs to Plot',roiList));
if (isempty(selectedROIs)), fprintf(1,'[%s]: No ROIs selected', mfilename); drawnow; return; end

% thresholds
vethresh = viewGet(vw,'cothresh');

for roi = selectedROIs
    roiCoords = vw.ROIs(roi).coords;
    [~, roiInd] = intersectCols(vw.coords, roiCoords);
    cVarExp = ve(ccm.rss(roiInd), ccm.rawrss(roiInd));
    rVarExp = ve(rm.rss(roiInd), rm.rawrss(roiInd));
    ii = cVarExp > vethresh & rVarExp > vethresh;
    ci = bootci(1000, {@mean,(cVarExp(ii)-rVarExp(ii))}, 'type', 'per');
    m = mean(ci);
    e = m-ci(1);
    fprintf(1,'[%s]: %0.2f +/- %0.2f %%\n', vw.ROIs(roi).name, m, e);
end


return

function varexp = ve(ssr,sst)
    varexp = 1 - (ssr ./ sst);
    varexp(~isfinite(varexp)) = 0;
    varexp = max(varexp, 0);
    varexp = min(varexp, 1);
return