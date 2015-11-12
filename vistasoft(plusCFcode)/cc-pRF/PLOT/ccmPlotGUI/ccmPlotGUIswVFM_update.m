function ccmPlotGUIswVFM_update(vw, M)
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
% M = ccmSWin(vw,M,voxel);

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

axes(M.ui.tsAxes); cla; hold on; %#ok<MAXES>

hTSeries = plot(M.x, M.tSeries(:,voxel), 'r-');
hFit = plot(M.x, pred(:, 1), 'b-', 'LineWidth', 1);
hResidual = plot(M.x, abs(M.tSeries(:,voxel)-pred(:,1)), 'c:');

allPlotted = [M.tSeries(:,voxel); pred(:,1); M.tSeries(:,voxel)-pred(:,1)];
axis([min(M.x) max(M.x) min(allPlotted) max(allPlotted)]);
h = axis;

for n=1:numel(M.sepx),
    plot([1 1].*M.sepx(n), [h(3) h(4)], 'y-');
end;

xlabel('Time (sec)','FontSize', 14);
ylabel('BOLD signal change (%)','FontSize', 14);

set(gca, 'FontSize', 14);
set(gca,'LineWidth',1);
set(gca,'color','black');
grid on
set(gca, 'XColor', [0.2 0.2 0.2]);
set(gca, 'YColor', [0.2 0.2 0.2]);
set(gcf,'color',[0.94 0.94 0.94]);

leg = legend('Voxel time series','Prediction','Residuals','Location','NorthWest');
set(leg, 'EdgeColor', 'black');
set(leg,'TextColor',[0.94 0.94 0.94]);
set(leg,'color','none');

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

% txt = sprintf('Variance explained VFM: %.2f%%', M.veR*100);
% set(M.ui.r2TextRef, 'String', txt);

% print out coords of voxel
txt = ['Target voxel coords: ' num2str(M.params.roi.coords(:,voxel)')];
set(M.ui.coordsText, 'String', txt);
M.target = num2str(M.params.roi.coords(:,voxel));

txt = ['Source CF coords: ' num2str([rfParams(1) rfParams(2) rfParams(3)])];
set(M.ui.xyzText, 'String', txt);

% txt = sprintf('Eccentricity RS = %.1f%s', M.et, ' deg');
% set(M.ui.eccText, 'String', txt);

% txt = sprintf('Eccentricity VFM = %.1f%s', M.eR, ' deg');
% set(M.ui.eccTextRef, 'String', txt);

% txt = sprintf('Polar angle RS = %.1f%s', M.pt, ' rad');
% set(M.ui.polText, 'String', txt);

% txt = sprintf('Polar angle VFM = %.1f%s', M.pR, ' rad');
% set(M.ui.polTextRef, 'String', txt);

% txt = sprintf('Sigma = %.1f%s', rfParams(4), ' mm');
% txt = sprintf('Sigma RS = %.1f%s', M.st, ' mm');
% set(M.ui.sigmaText, 'String', txt);

% txt = sprintf('Sigma VFM = %.1f%s', M.sR, ' mm');
% set(M.ui.sigmaTextRef, 'String', txt);


% from rmPlot: cache the data for the current time series plot as well
M.currTsData.x    = M.x;
M.currTsData.pred = pred;
M.currTsData.raw  = M.tSeries(:,voxel);
M.currTsData.sepx = M.sepx;
M.currTsData.pRF  = RF;

%% mod 4 nov 2014

axes(M.ui.cfAxes); cla; hold on; %#ok<MAXES>

% xAxis=240;
% yAxis=15;
% varex = plot(M.t,M.ve*10,'r-');
% ecc = plot(M.t,M.e,'c-');
% %stairs(mod(p+pi,2*pi),'g');
% pol = plot(M.t,M.p,'g-');
% sigma = plot(M.t,M.s,'b-');
% 
% % run = strcat('RS ',num2str(RS));
% % title(run,'FontSize', 14);
% xlim([0 xAxis]);
% ylim([0 yAxis]);
% xlabel('TR (1.5sec)','FontSize', 14);
% ylabel('relative units','FontSize', 14);
% 
% set(gca, 'FontSize', 14);
% set(gca,'LineWidth',1);
% set(gca,'color','black');
% grid on
% set(gca, 'XColor', [0.2 0.2 0.2]);
% set(gca, 'YColor', [0.2 0.2 0.2]);
% set(gcf,'color',[0.94 0.94 0.94]);
% 
% leg = legend('Corrected VE (x10)','Eccentricity (deg)','Polar angle (rad)','Sigma (mm)','Location','NorthWest');
% set(leg, 'EdgeColor', 'black');
% set(leg,'TextColor',[0.94 0.94 0.94]);
% set(leg,'color','none');



% store the updated M struct in the figure's user data
M.prevVoxel = voxel;  % also update the most recent voxel
set(M.fig, 'UserData', M);



if get(M.ui.showCheck, 'Value') == 1 && any(rfParams)
    ccmShowConnectiveFieldsw(vw, rfParams,M, M.distances,M.params.roi.coords(:,voxel));
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
