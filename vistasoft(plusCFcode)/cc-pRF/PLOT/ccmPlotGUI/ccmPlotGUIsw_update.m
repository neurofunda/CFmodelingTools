function ccmPlotGUIsw_update(vw, M)
% Refresh the cortico-cortical model plotting GUI.
%
% 2011 KVH: adapted from ccmPlotGUI_update
% 2014 NG: modified it

if ~exist('M','var') || isempty(M), M = get(gcf, 'UserData'); end

% get current voxel from the GUI
voxel = get(M.ui.voxel.sliderHandle, 'Value');

% get needed params
coords = M.coords(:,voxel);
if isequal(M.roi.viewType, 'Gray')  % convert coords into an index
    coords = M.coords(voxel);
end

M.modelNum = get(M.ui.model, 'Value');

%% load sliding windows data
M = ccmSWin(M,voxel);

% compute RF, get tSeries for this voxel
[pred RF rfParams variance_explained CFtSeries] = ccmPlotGUI_makePrediction(M, coords, voxel);

% check if the prediction for this voxel is empty. If so, don't plot it:

if all( isnan(pred) | isinf(pred) )
    axes(M.ui.tsAxes);  cla;   %#ok<MAXES>
    AX = axis;
    text(AX(2) + .5*(AX(2)-AX(1)), AX(3) + .5*(AX(4)-AX(3)), ...
        '(No data available for this voxel)', ...
        'FontSize', 12, ...
        'HorizontalAlignment', 'center');
    return
end

% store the current prediction
M.prediction = pred;
M.RF = RF;

% roiList = viewGet(vw, 'roinames');
% inputFolderWholeVFM = '/Users/Nicolas/Desktop/subj_3/WindowedRS/Analysis/VFM/';
% dataVFM = strcat(inputFolderWholeVFM,'cfData_VisualAverages_',roiList{1},'_',roiList{2},'.mat');
% load(dataVFM);
% Source
[c, roiInd] = intersectCols(vw.coords, vw.ROIs(1).coords);
tSeriesSource  = loadtSeries(vw, 1);
tSeriesSource = tSeriesSource(:,roiInd(:));
tSeriesSource = lowpass(tSeriesSource);

axes(M.ui.tsAxes); cla; hold on; %#ok<MAXES>


hTSeries = plot(M.x, M.tSeries(:,voxel), 'r-');
hFit = plot(M.x, pred(:, 1), 'k--', 'LineWidth', 1);
hResidual = plot(M.x, (M.tSeries(:,voxel)-pred(:,1)), 'k:');
CFsignal = plot(M.x, CFtSeries(:,M.modelNum), 'b-');

allPlotted = [M.tSeries(:,voxel); pred(:,1); M.tSeries(:,voxel)-pred(:,1);CFtSeries(:,M.modelNum)];
axis([min(M.x) max(M.x) min(allPlotted) max(allPlotted)]);
h = axis;

% for n=1:numel(M.sepx),
%     plot([1 1].*M.sepx(n), [h(3) h(4)], 'y-');
% end;

xlabel('Time (sec)','FontSize', 12);
ylabel('BOLD signal change (%)','FontSize', 12);

set(gca, 'FontSize', 12);
set(gca,'LineWidth',1);
%set(gca,'color','black');
%grid on
set(gca, 'XColor', [0.2 0.2 0.2]);
set(gca, 'YColor', [0.2 0.2 0.2]);
set(gcf,'color',[0.94 0.94 0.94]);

leg = legend('Voxel time series','Prediction','Residuals','CF time series','Location','NorthWest');
%set(leg, 'EdgeColor', 'black');
legend(leg,'boxoff')
set(leg,'TextColor',[0 0 0]);
set(leg,'color','none');


% set the user data of each checkbox to point to the PLOT curves
% we just created. The callbacks for these checkboxes will then toggle
% the visibility of each curve.
set(M.ui.tsCheck, 'UserData', hTSeries);
set(M.ui.predCheck, 'UserData', hFit);
set(M.ui.resCheck, 'UserData', hResidual);
set(M.ui.showCFsignal, 'UserData', CFsignal);
if get(M.ui.tsCheck, 'Value')==0, set(hTSeries, 'Visible', 'off'); end
if get(M.ui.predCheck, 'Value')==0, set(hFit, 'Visible', 'off'); end
if get(M.ui.resCheck, 'Value')==0, set(hResidual, 'Visible', 'off'); end
if get(M.ui.showCFsignal, 'Value')==0, set(CFsignal, 'Visible', 'off'); end


% (2) set text fields

% print out percent variance (R^2) explained
% txt = sprintf('Variance explained RS: %.2f%%', variance_explained*100);
txt = sprintf('Variance explained RS: %.2f%%',  M.vetRS*100);
set(M.ui.r2Text, 'String', txt);

txt = sprintf('Variance explained VFM: %.2f%%', M.vetVFM*100);
set(M.ui.r2TextRef, 'String', txt);

% print out coords of voxel
txt = ['Target voxel coords: ' num2str(M.params.roi.coords(:,voxel)')];
set(M.ui.coordsText, 'String', txt);
M.target = num2str(M.params.roi.coords(:,voxel));

txt = ['CF center coords: ' num2str([rfParams(1) rfParams(2) rfParams(3)])];
set(M.ui.xyzText, 'String', txt);

txt = sprintf('Eccentricity RS = %.1f%s', M.etRS, ' deg');
set(M.ui.eccText, 'String', txt);

txt = sprintf('Eccentricity VFM = %.1f%s', M.etVFM, ' deg');
set(M.ui.eccTextRef, 'String', txt);

txt = sprintf('Polar angle RS = %.1f%s', M.ptRS, ' rad');
set(M.ui.polText, 'String', txt);

txt = sprintf('Polar angle VFM = %.1f%s', M.ptVFM, ' rad');
set(M.ui.polTextRef, 'String', txt);

txt = sprintf('Sigma RS = %.1f%s', M.stRS, ' mm');
set(M.ui.sigmaText, 'String', txt);

txt = sprintf('Sigma VFM = %.1f%s', M.stVFM, ' mm');
set(M.ui.sigmaTextRef, 'String', txt);


% from rmPlot: cache the data for the current time series plot as well
M.currTsData.x    = M.x;
M.currTsData.pred = pred;
M.currTsData.raw  = M.tSeries(:,voxel);
M.currTsData.sepx = M.sepx;
M.currTsData.pRF  = RF;

%% mod 4 nov 2014
axes(M.ui.cfAxes); cla; hold on; %#ok<MAXES>

xAxis=240;
yAxis=15;
varex = plot(M.t,M.ve*10,'r-');
ecc = plot(M.t,M.e,'b-');
pol = plot(M.t,M.p,'m-');
sigma = plot(M.t,M.s,'k-');


% run = strcat('RS ',num2str(RS));
% title(run,'FontSize', 12);
xlim([0 xAxis]);
ylim([0 yAxis]);
xlabel('TR (1.5sec)','FontSize', 12);
ylabel('relative units','FontSize', 12);

set(gca, 'FontSize', 12);
set(gca,'LineWidth',1);
%set(gca,'color','black');
%grid on
set(gca, 'XColor', [0.2 0.2 0.2]);
set(gca, 'YColor', [0.2 0.2 0.2]);

set(gcf,'color',[0.94 0.94 0.94]);

leg = legend('Corrected VE (x10)','Eccentricity (deg)','Polar angle (rad)','Sigma (mm)','Location','NorthWest');
%set(leg, 'EdgeColor', 'black');
legend(leg,'boxoff')
set(leg,'TextColor',[0 0 0]);
set(leg,'color','none');

%% mod 26 march 2015
% set WTC plot
 axes(M.ui.corrAxes); cla; hold on; %#ok<MAXES>
% if get(M.ui.showCFsignal, 'Value')==1, wtc(M.tSeries(:,voxel),CFtSeries(:, 1),'Dj',1/6,'mcc',0,'ad',[20 20],'as',0.5,'ahs',0.2);
% else wtc(M.tSeries(:,voxel),pred(:, 1),'Dj',1/6,'mcc',0,'ad',[20 20],'as',0.5,'ahs',0.2); end
plot(M.t,M.cfCorr_RS,'b');
plot(M.t,M.cfCorr_VFM,'r');
xlim([0 xAxis]);
ylabel('Correlation coefficient','FontSize', 12);
leg = legend('Correlation with RS-CF center','Correlation with VFM-CF reference center','Location','SouthWest');
legend(leg,'boxoff')
set(gca, 'FontSize', 12);
set(leg,'TextColor',[0 0 0]);
set(leg,'color','none');





% store the updated M struct in the figure's user data
M.prevVoxel = voxel;  % also update the most recent voxel

set(M.fig, 'UserData', M);


if get(M.ui.showConnections, 'Value') == 1 && any(rfParams)
    
    %% Another way too look at the CF "connections" is to extract ROI referred indices and plot lines between source and target indices, also the size is visualized as circles and the corrected EV as color.
    % gather targetROI referred sampling size
    targetROI = viewGet(vw, 'CurrentROI');
    targetCoords = vw.ROIs(2).coords;
    %% obtain ROI indices
    [c, roiInd] = intersectCols(vw.coords, targetCoords);
    %% define CF model structure
    ccm = vw.ccm.models{vw.ccm.modelNum};
    %% get model indices
    ind = ccmGet(ccm,'indices');
    %% get ROI referred marameters for the CF model
    sigma = ccm.sigma.major(ind);
    rss = ccm.rss(ind);
    rawrss = ccm.rawrss(ind);
    npoints = ccm.npoints;
    ntrends = ccm.ntrends; %% you must set this to 0 if not detrending
    %ntrends = 0;
    ve = 1 - (rss ./ rawrss); % variance explaned
    cVE = ve-(1-ve)*(ntrends/(npoints-ntrends-1)); % corrected variance explained
    %% get CF model coordinates
    x0 = ccmGet(ccm,'x0'); y0 = ccmGet(ccm,'y0'); z0 = ccmGet(ccm,'z0');
    centerCoords = cat(2,x0(ind)',y0(ind)',z0(ind)')';
    sourceCoords = vw.ROIs(1).coords;
    %% transform coordinates into ROI indices,
    for i=1:length(centerCoords)
        [cSource(:,i), indxSource(i)] = intersectCols(sourceCoords, centerCoords(:,i));
        con(i) = indxSource(i);
    end
    
    %% Load reference VFM CF connections
    roiList = viewGet(vw, 'roinames');
    inputFolderWholeVFM = '/Users/Nicolas/Desktop/subj_3/WindowedRS/Analysis/VFM/';
    dataVFM = strcat(inputFolderWholeVFM,'cfData_VisualAverages_',roiList{1},'_',roiList{2},'.mat');
    load(dataVFM);
    
    %% Plot connections between indices, highlight current and reference connection
    connections = figure(69);
    clf(connections);
    tIndex = [1:length(centerCoords)];
    CFindex = zeros(length(tIndex));
    %scatter(CFindex(1,:), tIndex,50,((veVFM-min(veVFM(:)))./(max(veVFM(:)-min(veVFM(:))))));
    scatter(CFindex(1,:), tIndex,50);
    colormap(jet);
    caxis([0 2*pi]);
    hold on
    cmap = colormap(jet(length(cVE)));
    for i=1:length(tIndex)
        patch([1,0],[con(i),tIndex(i)],'k', 'EdgeAlpha', 0.1, 'FaceColor', 'none'); % RS CF connections
        %patch([1,0],[cf.conIndex(i),tIndex(i)],'k', 'EdgeAlpha', 0.1, 'FaceColor', 'none'); % VFM CF reference connections
    end
    patch([1,0],[con(voxel),tIndex(voxel)],'r','LineWidth',2,'EdgeColor', 'r'); % RS current voxel connection
    patch([1,0],[cf.conIndex(voxel),tIndex(voxel)],'b','LineWidth',2,'EdgeColor', 'b'); % VFM current voxel connection
    targetVoxel = ones(length(con));
    S = sigma*15*pi.^2;
    scatter(targetVoxel(1,:), con,S,cVE);
    colormap(jet);
    caxis([0 1])
    h = colorbar;
    title(h,'Variance explained')
    xlim([0 1.2]);
    set(gcf, 'color', 'w');
    %set(gcf,'units','centimeters','position',[0 0 400 800])
    set(gca,'LineWidth',1)
    set(gca, 'FontSize', 14);
    set(gca,'xtick',[0 1 ])
    set(gca,'xticklabel',{'Target #','Source #'})
    ylabel('Voxel ');
    title('CF connections ');
    
    
    %% Load and preprocess data
    % Source
    [c, roiInd] = intersectCols(vw.coords, vw.ROIs(1).coords);
    tSeriesSource  = loadtSeries(vw, 1);
    tSeriesSource = tSeriesSource(:,roiInd(:));
    % Convert to percent bold
    tSeriesSource  = raw2pc(tSeriesSource(:,cf.conIndex(voxel)));
    % Low pass filter
    tSeriesSource = lowpass(tSeriesSource);
    % Detrend
    trends  = rmMakeTrends(M.params);
    b = pinv(trends)*tSeriesSource;
    tSeriesSource = tSeriesSource - trends*b;
    % Target
    [c, roiInd] = intersectCols(vw.coords, vw.ROIs(2).coords);
    tSeriesTarget  = loadtSeries(vw, 1);
    tSeriesTarget = tSeriesTarget(:,roiInd(:));
    % Convert to percent bold
    tSeriesTarget  = raw2pc(tSeriesTarget(:,voxel));
    % Low pass filter
    tSeriesTarget = lowpass(tSeriesTarget);
    % Detrend
    trends  = rmMakeTrends(M.params);
    b = pinv(trends)*tSeriesTarget;
    tSeriesTarget = tSeriesTarget - trends*b;
    
    %     histo = figure (1000)
    %     clf(histo);
    %     hist(M.tSeries(:,voxel));
    %     set(gcf, 'color', 'w');
    
    
    %% Time Series Analysis Figure Current connection
    currentCF = figure(100);
    clf(currentCF);
    %title('Time Series analysis for current RS CF connection');
    
    subplot ('Position',[0.8 0.75  0.1 0.2]);  
    
    [Cxy,F] = mscohere(CFtSeries(:, M.modelNum),M.tSeries(:,voxel),hamming(100),75,[],1/1.5);
%     df = 2 * length(M.tSeries(:,voxel))/80; %degrees of freedom
%     error_coherence = sqrt(2)*(1-Cxy)./(abs(sqrt(Cxy))*sqrt(df)); % equation 9.82
%     
    %[pks,locs] = findpeaks(double(Cxy),'MinPeakHeight',0.75);
    [pks,locs] = max(Cxy);
    plot(F,Cxy,'-k');
    hold
    plot(F(locs),pks,'p','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',20);
    xlim([ 0 0.1]); xlabel('Frequency (Hz)'); 
    ylim([ 0 1]); ylabel('Magnitude');
    title('Magnitude-Squared Coherence');
    %grid on    
    view(90,-90);
    
    colormap('winter');
    subplot('Position',[0.13 0.75  0.58 0.2]);
    %axes('Position',[0.12 0.71  0.6 0.25])
    %wtc(M.tSeries(:,voxel),CFtSeries(:, 1),'Dj',1/10,'mcc',0,'ad',[20 20],'as',0.5,'ahs',0.2);
    %title('Cross-Wavelet Coherence for current RS CF connection');
    [C,phi,S12,S1,S2,t,f,confC,phistd] = cohgramc(M.tSeries(:,voxel),CFtSeries(:, M.modelNum),[100,25],...
        struct('tapers', [2 3],'pad', -1,'fpass', [0 0.1],'err', [1,0.05],'Fs',[1/1.5]));
    
    fInterest = F(locs);    
    fbin = find(f>fInterest-0.005 & f<fInterest+0.005); %% bin of interest
    
    C(C<confC)=NaN;
    %  Note that phi + 2 phistd and phi - 2 phistd will give 95% confidence
    phi(phistd>.5)=NaN;
    pcolor(t,f, C'); shading flat;
    hold on
    phaseplot(t,f(:,fbin),phi(:,fbin),1/25,1/25);
    coh = colorbar; %set(get(coh,'title'),'String','Coherence');
    caxis([0 1]);
    title('Spectral Coherence with RS-derived CF center');
    xlabel('Time (sec)'); ylabel('Frequency (Hz)');
    
   
    
    subplot(4,2,[3,4]);
    plot(M.x,M.tSeries(:,voxel),'r');
    hold on;
    plot(M.x,CFtSeries(:, M.modelNum),'b');
    plot(M.x, pred(:, 1),'--k');
    %grid on
    %title('Time Series');
    xlim([ 0 360]);
    xlabel('Time(s)');
    ylabel('BOLD signal change (%)')
    legend('Target voxel','CF center','CF prediction');
    leg = legend(leg,'boxoff')
    set(gca, 'FontSize', 12);
    set(leg,'TextColor',[0 0 0]);
    set(leg,'color','none');
    
    subplot(4,2,[5,6]);
    plot(M.t,M.cfCorr_RS,'k');
    grid on
    title('Windowed Correlation with RS-derived CF center');
    xlim([0 240]);
    xlabel('TR: 1.5 sec');
    ylabel('Correlation coefficient');
    ylim([-0.25 1]);
%     leg = legend('Correlation with RS-derived CF center','Location','SouthEast');
%     legend(leg,'boxoff')
%     set(gca, 'FontSize', 12);
%     set(leg,'TextColor',[0 0 0]);
%     set(leg,'color','none');
    
    subplot(4,2,[7,8])
    %figure,
    [xycov,lags] = xcov(M.tSeries(:,1),CFtSeries(:, M.modelNum));
    %for i=1:616
    %[xycov,lags] = xcov(M.tSeries(:,i),CFtSeries(:,i));
    plot(lags,xycov/max(xycov),'-k');
    %hold on
    %end
    %grid on
    ylim([-0.5 1]);
    title('Cross-Correlation with RS-derived CF center');
    xlabel('Lag (TR: 1.5sec)');
    ylabel('Magnitude');
    set(gcf, 'color', 'w');
    
    
    
    
%     %% Time Series Analysis Figure Reference connection
%     referenceCF = figure(200);
%     clf(referenceCF);
%     %title('Time Series analysis for  reference VFM CF connection on RS data');
%     
%     
%     subplot ('Position',[0.8 0.75  0.1 0.2]); 
%     
%     [Cxy,F] = mscohere(tSeriesTarget,tSeriesSource,hamming(100),75,[],1/1.5);
%     %[pks,locs] = findpeaks(double(Cxy),'MinPeakHeight',0.75);
%     [pks, locs] = max(Cxy);
%     plot(F,Cxy,'-k');
%     hold
%     plot(F(locs),pks,'p','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',20);
%     xlim([ 0 0.1]); xlabel('Frequency (Hz)');
%     ylim([ 0 1]); ylabel('Magnitude');
%     title('Magnitude-Squared Coherence');
%     %grid on
%     view(90,-90);
%     set(gcf, 'color', 'w');
%     
%   
%     colormap('winter');
%     subplot('Position',[0.13 0.75  0.58 0.2]);
%     %wtc(tSeriesTarget,tSeriesSource,'Dj',1/10,'mcc',0,'ad',[20 20],'as',0.5,'ahs',0.2);
%     %title('Cross-Wavelet Coherence for  reference VFM CF connection on RS data');
%     %xlabel('TR: 1.5 sec');
%     %ylabel('Period (sec)');
%     [C,phi,S12,S1,S2,t,f,confC,phistd] = cohgramc(tSeriesTarget,tSeriesSource,[100,25],...
%         struct('tapers', [2 3],'pad', -1,'fpass', [0 0.1],'err', [1,0.05],'Fs',[1/1.5]));
%     
%     
%     fInterest = F(locs);    
%     fbin = find(f>fInterest-0.005 & f<fInterest+0.005); %% bin of interest
%     
%     
%     C(C<confC)=NaN;
%     %  Note that phi + 2 phistd and phi - 2 phistd will give 95% confidence
%     phi(phistd>.5)=NaN;
%     pcolor(t,f, C'); shading flat;
%     hold on
%     phaseplot(t,f(:,fbin),phi(:,fbin),1/25,1/25);
%     coh = colorbar; %set(get(coh,'title'),'String','Coherence');
%     caxis([0 1]);
%     title('Cross-Coherence with VFM-derived CF center');
%     xlabel('Time (sec)'); ylabel('Frequency (Hz)');
% 
%     
%     subplot(4,2,[3,4]);
%     plot(tSeriesTarget,'r');
%     hold on;
%     plot(tSeriesSource,'b');
%     %grid on
%     %title('Time Series');
%     xlabel('TR: 1.5 sec');
%     ylabel('BOLD signal change (%)');    
%     legend('Target voxel','CF center');
%     leg = legend(leg,'boxoff')
%     set(gca, 'FontSize', 12);
%     set(leg,'TextColor',[0 0 0]);
%     set(leg,'color','none');
% 
%     
%     subplot(4,2,[5,6]);
%     plot(M.t,M.cfCorr_VFM,'k');
%     grid on
%     title('Windowed Correlation with VFM-derived CF center');
%     xlabel('TR: 1.5 sec');
%     ylabel('Correlation coefficient');
%     ylim([-0.5 1]);
%     %leg = legend('Windowed correlation with VFM-derived CF center','Location','SouthEast');
%     %legend(leg,'boxoff')
%     %set(gca, 'FontSize', 12);
%     %set(leg,'TextColor',[0 0 0]);
%     %set(leg,'color','none');
%     
%     
%     subplot(4,2,[7,8]);
%     [xycov,lags] = xcov(tSeriesTarget,tSeriesSource);
%     plot(lags,xycov/max(xycov),'-k');
%     %grid on
%     ylim([-0.5 1]);
%     title('Cross-Covariance with VFM-derived CF center');
%     xlabel('Lag (TR: 1.5 sec)');
%     ylabel('Magnitude');
%     set(gcf, 'color', 'w');
    

    clear cf
    
    
end

if get(M.ui.showCheck, 'Value') == 1 && any(rfParams)
    ccmShowConnectiveFieldsw(vw, rfParams,M, M.distances,M.params.roi.coords(:,voxel),centerCoords(:,voxel));
end


if get(M.ui.showMov, 'Value') == 1 && any(rfParams)
    %% video
    name = sprintf('Polar Map %s (%s)', M.roi.name, M.dataType);
    M.polMap = figure('Color','white', 'Name', name, ...
        'Units', 'pixels', 'Position', [50 0 300 300]);
    visualFieldCFmov(M,voxel);
    % set(M.polMap,'UserData', M);
end

return;

function tSeries = lowpass(tSeries)
TR = 1.5;
fc = 0.1;
fs = 1/TR;
[b, a] = butter(4,fc/(fs/2),'low');
tSeries = filtfilt(b,a,double(tSeries));
return

%---------------------------------
function data=raw2pc(data)
dc   = ones(size(data,1),1)*mean(data);
data = ((data./dc) - 1) .*100;
return;
%---------------------------------




