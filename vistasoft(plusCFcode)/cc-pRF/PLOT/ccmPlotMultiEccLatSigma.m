function ccmPlotMultiEccLatSigma(vw, ROIlist)
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

% check ROIs
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

prompt = {'Enter binsize:'};
dlg_title = 'Please provide binsize';
num_lines = 1;
def = {'0.25','hsv'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
binsize = eval(answer{1});

prompt = {'Enter eccentricity range:'};
dlg_title = 'Please provide eccentricity range';
num_lines = 1;
def = {['[0 ' num2str(max([rm.params.stim(:).stimSize])) ']'],'hsv'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
range = eval(answer{1});

% set up plot
graphwin = figure; hold on;
set(graphwin, 'Color', 'w')
c = jet(nROIs)*.9;

% initialize a legend
legendtxt = cell(1,nROIs);

% initialize data struct
data = cell(1, nROIs); 

% suppress individual plots from calls to rmPlotEccSigma
plotFlag = false; 

% loop through ROIs
for ii = 1:nROIs
    vw = viewSet(vw, 'curroi', selectedROIs(ii));
    data{ii} = ccmPlotEccLatSigma(vw, rm, plotFlag, binsize, range);
    data{ii}.roi = viewGet(vw, 'roiname');
    legendtxt{ii} = data{ii}.roi;
    figure(graphwin);
    plot(data{ii}.xfit, data{ii}.yfit, '-', 'color', c(ii,:), 'LineWidth', 2)
end

legend(legendtxt,2);
maxy = 0;

% add the data points for each plot
for ii = 1:nROIs
    if ~exist('bootstrp','file')
        errorbar(data{ii}.x,data{ii}.y,data{ii}.ysterr, 'x', 'color', c(ii,:));   
    else
        p1 = plot(data{ii}.b_xfit,data{ii}.b_upper,'--','LineWidth',1.5);
        p2 = plot(data{ii}.b_xfit,data{ii}.b_lower,'--','LineWidth',1.5);
        set(p1,'Color',[.6 .6 .6]);
        set(p2,'Color',[.6 .6 .6]);
        plot(data{ii}.x,data{ii}.y,'ko',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',c(ii,:),...
            'MarkerSize', 8);
        maxy = max(maxy,max(data{ii}.y));
    end
end
axis([range(1) range(2) 0 maxy+1]);
ylabel('Connective Field Radius (mm)');
xlabel('Eccentricity (deg)');
title(['Sampling from ' vw.ccm.params.sourceROI]);

return
