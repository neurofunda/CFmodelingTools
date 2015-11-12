function ccmPlotMultiSigma(vw, ROIlist)
% Plots the mean cc-pRF size within a specified range of eccentricities
%
% 2010 KVH: wrote it

% check view struct
if notDefined('vw'), vw = getCurView; end

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

% select ROIs 
if (notDefined('ROIlist'))
    roiList=viewGet(vw, 'roinames');
    if isempty(roiList)
        fprintf(1,'[%s]: please load and select at least one ROI for plotting\n', ...
        mfilename); drawnow; return;
    end
    selectedROIs = find(buttondlg('ROIs to Plot',roiList));
elseif ROIlist == 0,
    selectedROIs = 1:length(viewGet(vw, 'ROIs'));
else
    selectedROIs=ROIlist;
end

nROIs=length(selectedROIs);
if (nROIs==0), fprintf(1,'[%s]: No ROIs selected', mfilename); drawnow; return; end

% Load the retinotopic model file
rmFile = getPathStrDialog(dataDir(vw),...
        'Choose pRF model file name', 'retModel*.mat');    
rm = load(rmFile);

% specify the eccentricity range
prompt = {'Enter eccentricity range:'};
dlg_title = 'Please provide eccentricity range';
num_lines = 1;
def = {['[0 ' num2str(max([rm.params.stim(:).stimSize])) ']'],'hsv'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
range = eval(answer{1});

% set up plot
graphwin = figure; hold on;
set(graphwin, 'Color', 'w')

% initialize
m = zeros(1, nROIs); 
e = zeros(1, nROIs);
ROInames = cell(1,nROIs);

% suppress individual plots from calls to ccmPlotEccSigma
plotFlag = false; 

% loop through ROIs
fprintf(1,'[%s]: Computing 95%% bootstrap confidence intervals... \n', mfilename);
drawnow;
for ii = 1:nROIs
    vw = viewSet(vw, 'curroi', selectedROIs(ii));
    data = ccmPlotEccSigma(vw, rm, plotFlag, .1, range);
    wmean = @(x,w) sum(x.*w)./sum(w);
    ci = bootci(1000,{wmean,data.sig,data.ve},'type','per');
    m(ii) = mean(ci);
    e(ii) = m(ii) - ci(1);
    ROInames{ii} = vw.ROIs(selectedROIs(ii)).name; 
    fprintf(1,'[%s]: %s: sigma =%6.3f +/-%6.3f \n', mfilename, ROInames{ii}, m(ii), e(ii));
    drawnow;
end

bh = bar(m);
ch=get(bh,'children');
cd=repmat(1:numel(m),5,1);
cd=[cd(:);nan];
set(ch,'facevertexcdata',cd);
colormap(jet(128));
errorbar(m,e,'k','LineStyle','none');
ylabel('Connective Field Radius (mm)'); ylim([0 15])
set(gca,'XTick',1:nROIs,'XTickLabel',ROInames);
title(['Sampling from ' vw.ccm.params.sourceROI]);

return
