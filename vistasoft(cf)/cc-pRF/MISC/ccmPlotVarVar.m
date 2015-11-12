function difference = ccmPlotVarVar(vw, rm, pf)
% cortico-cortical model must be loaded
mrGlobals;

try
	ccm = vw.ccm.models{vw.ccm.modelNum};
catch %#ok<CTCH>
    ccm = [];
end;

if isempty(ccm)
    fprintf(1,'[%s]: please load a cortico-cortical model first\n', ...
        mfilename); drawnow; return;
end

roi = vw.selectedROI;
if isequal(roi,0),
    fprintf(1,'[%s]: please load and select a ROI for plotting\n', ...
        mfilename); drawnow; return;
end

if notDefined('rm')
    rmFile = getPathStrDialog(dataDir(vw),...
        'Choose pRF model file name', 'retModel*.mat');    
    rm = load(rmFile);
end

if notDefined('pf'), pf = 1; end

targetROI = viewGet(vw, 'CurrentROI');
 
targetCoords = vw.ROIs(targetROI).coords;
[~, roiInd] = intersectCols(vw.coords, targetCoords);

ssr = ccm.rss(roiInd);
sst = ccm.rawrss(roiInd);
ccm.varexp = varexp(ssr,sst);

ssr = rm.model{1}.rss(roiInd);
sst = rm.model{1}.rawrss(roiInd);
rm.varexp = varexp(ssr,sst);

% thresholds
thresh.varexp = 0;

% find useful data given thresholds
ii = ccm.varexp > thresh.varexp & rm.varexp > thresh.varexp; 

% get the difference
difference = ccm.varexp(ii) - rm.varexp(ii); 
                                
% plot if requested
if pf
    figure('Color', 'w'); hold on;
    plot(rm.varexp(ii),ccm.varexp(ii),'.','Color',[.4 .4 .4])
    plot([0 1], [0 1],'Color','r', 'LineWidth',2, 'LineStyle', '--');
    axis([0 1 0 1])
    axis square
else
    save(['ccmVarExp_' vw.ROIs(targetROI).name '-' dataTYPES(viewGet(vw,'curdt')).name '.mat'], 'difference');
end



return

function ve = varexp(ssr,sst)
    ve = 1 - (ssr ./ sst);
    ve(~isfinite(ve)) = 0;
    ve = max(ve, 0);
    ve = min(ve, 1);
return
