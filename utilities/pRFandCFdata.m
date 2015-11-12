function [ ] = pRFandCFdata(view,firstDatatype,lastDatatype,VEthr,pRFmodel,nFrames,TR)
%--------------------------------------------------------------------------
% This function retrieves already computed pRF model parameters, the times series of a
% source and a target ROI (with the option to filter) and compute CF
% models. 
% 
% The function depends on the following inputs:
% 
% -     A loaded mrVista session (view) with loaded source and target ROIs
% -     A folder "Analysis" inside the mrVista session folder
% (mkdir('Analysis'). The results will be saved there.
% -     Path to a folder containing the pRF model
% -     Datatype numbers (This is vestigial for now. It can be useful to compute
%       many many many dataTypes at once.)
% -     Variance explained (threshold for plotting pRF maps)
% -     TR

% The output of this function is the structure pRFandCFdata_dataType_sourceROI_targetROI, which contains:
%
% -   Source and target ROI time series
% -   Associated pRF variance explained (VE)
% -   Associated pRF size (Sigma)
% -   Associated Visual field coordinates (Eccentricity and Polar angle)
% -   Associated Cortical coordinates
% -   Cortical distances between all pair of nodes (voxels)
% -   The structure"cf", which contains the CF parameters:
%   -   Variance explained
%   -   Corrected variance explained
%   -   Sigma
%   -   Rss, Raw Rss, beta, npoints, ntrends
%   -   CF center coordinates
%   -   Cartesian and polar visuotopic coordinates
%   -   Target voxel coordinates
%   -   Connection index (each target ROI index and its associated CF index in the source ROI)
% 
% This function also computes the power spectrum of the detrended time series, plots them and save the images as pdf files. 
% The power spectrum is computed using Chronux.
% 
% Example:
% 
% %% Analysis
% %% Load data into mrVista  (run the first part of the script in Matlab 2011b)
% % Clean memory
% close all
% clear all
% % Open mrVista
% mrVista 3
% view = VOLUME{1};
% mrGlobals;
% 
% define folders
% %folder = '/Users/visionlab/Dropbox (Vision)/RET-3T/subj_1/'; % Office
% folder='/Volumes/Data/Dropbox (Vision)/RET-3T/subj_1/'; % Home
% pRFfolder = [folder 'Gray/Averages/'];
% pRFmodel = strcat(pRFfolder,'retModel-20140723-130624-fFit','.mat');
% 
% %% Collect time series and coordinates for left hemisphere
% % Left hemisphere and ROIs
% hemis = 'L';
% source = 'V1';
% target = 'V3';
% % Load ROIS
% ROIS{1} = strcat(hemis,source) % Source ROI
% ROIS{2} = strcat(hemis,target) % Target ROI
% fprintf('Loading source ROI:%s\n',ROIS{1});
% fprintf('Loading target ROI:%s\n',ROIS{2});
% view = loadROI(view, ROIS,[], [], 0, 1);
% roiList = viewGet(view, 'roinames');
% 
% % If plotting 3D maps:
% % % Load mesh in workspace
% % mesh = strcat([folder 'Anatomy/'],'mesh(inflated)_', hemis,'.mat');
% % % Set view mapMode
% % view.ui.mapMode = setColormap(view.ui.mapMode, 'hsvTbCmap');
% % % Load mesh in mrMesh (open 3D window)
% % view = meshLoad(view, mesh, 1);
% % % Load mesh setting
% % msh = viewGet(view, 'selectedMesh');
% % meshRetrieveSettings(msh, 'lMesh');
% % Compute CF models
% 
% % VFM
% pRFandCFdata(view,4,4,0,pRFmodel,[],128,1.5);
%
% nicolas.gravel@gmail.com   11-11-2015 University of Groningen, The Netherlands.
%--------------------------------------------------------------------------

% The following parameters are needed for detrending with the function "rmMakeTrends"
params = struct('stim',struct('nFrames', nFrames,'framePeriod',TR,'nDCT',3,'nUniqueRep',1));
% Retrieve ROIs
roiList = viewGet(view, 'roinames');
% Load pRF model
rm = load(pRFmodel);

%--------------------------------------------------------------------------
% Source pRF coords
%--------------------------------------------------------------------------
sourceROI = roiList{1};
% ROI brain coordinates (3D)
sourceCoords = view.ROIs(1,1).coords;
[c, roiInd] = intersectCols(view.coords, sourceCoords);
x0_sourceROI = []; y0_sourceROI = [];
for j = 1:length(sourceCoords)
    x0_sourceROI = [x0_sourceROI rm.model{1}.x0(roiInd(j))];
    y0_sourceROI = [y0_sourceROI rm.model{1}.y0(roiInd(j))];
end
% convert to polar coordinates
[pol, ecc] = cart2pol(x0_sourceROI, y0_sourceROI);
pol = mod(pol, 2*pi);
% must do some sanity checks here
pol(pol == Inf)  = max(pol(isfinite(pol(:))));
pol(pol == -Inf) = min(pol(isfinite(pol(:))));
pol(isnan(pol)) = 0;
pol = max(pol,0);
pol = min(pol,2*pi);
ecc(ecc == Inf)  = max(ecc(isfinite(ecc(:))));
ecc(ecc == -Inf) = min(ecc(isfinite(ecc(:))));
ecc(isnan(ecc)) = 0;
% Visual field coordinates
sourceVCoords = cat(2,ecc',pol')';
% Variance explained
ve = 1 - (rm.model{1}.rss ./ rm.model{1}.rawrss);
ve(~isfinite(ve)) = 0;
ve = max(ve, 0);
ve = min(ve, 1);
sourceVE = ve(roiInd);
sourceSigma = rm.model{1}.sigma.major(roiInd);

%--------------------------------------------------------------------------
% target pRF coords
%--------------------------------------------------------------------------
targetROI = roiList{2};
% ROI brain coordinates (3D)
targetCoords = view.ROIs(1,2).coords; 
[c, roiInd] = intersectCols(view.coords, targetCoords);
x0_targetROI = []; y0_targetROI = [];
for j = 1:length(targetCoords)
    x0_targetROI = [x0_targetROI rm.model{1}.x0(roiInd(j))];
    y0_targetROI = [y0_targetROI rm.model{1}.y0(roiInd(j))];
end
% convert to polar coordinates
[pol, ecc] = cart2pol(x0_targetROI, y0_targetROI);
pol = mod(pol, 2*pi);
% must do some sanity checks here
pol(pol == Inf)  = max(pol(isfinite(pol(:))));
pol(pol == -Inf) = min(pol(isfinite(pol(:))));
pol(isnan(pol)) = 0;
pol = max(pol,0);
pol = min(pol,2*pi);
ecc(ecc == Inf)  = max(ecc(isfinite(ecc(:))));
ecc(ecc == -Inf) = min(ecc(isfinite(ecc(:))));
ecc(isnan(ecc)) = 0;
% Visual field coordinates
targetVCoords = cat(2,ecc',pol')';
% Variance explained
ve = 1 - (rm.model{1}.rss ./ rm.model{1}.rawrss);
ve(~isfinite(ve)) = 0;
ve = max(ve, 0);
ve = min(ve, 1);
targetVE = ve(roiInd);
targetSigma = rm.model{1}.sigma.major(roiInd);

%--------------------------------------------------------------------------
% Get distance matrix for source and target ROI
%--------------------------------------------------------------------------
rois = viewGet(view, 'rois');
sourceROI = rois(1);
targetROI = rois(2);
sourceDistances = ccmVoxToVoxDist(sourceROI, view, view.mmPerVox);
targetDistances = ccmVoxToVoxDist(targetROI, view, view.mmPerVox);

%--------------------------------------------------------------------------
% Loop trough datatypes
%--------------------------------------------------------------------------
global dataTYPES;
d = 0;
for j = firstDatatype:lastDatatype
    % Load datatypes
    d = d+1;
    
    view = selectDataType(view,j);
    curDataType = viewGet(view,'curDataType');
    fprintf('Loading Datatype:%s\n',dataTYPES(1,curDataType).name);
        
    %% Load source time series
    [c, roiInd] = intersectCols(view.coords, view.ROIs(1).coords);
    tSeriesSource  = loadtSeries(view, 1);
    tSeriesSource = tSeriesSource(:,roiInd(:,:));
    tSeriesSource= raw2pc(tSeriesSource);
    
    %--------------------------------------------------------------------------
    % Detrend
    %--------------------------------------------------------------------------
    trends  = rmMakeTrends(params);
    trends = single(trends);
    b = pinv(trends)*tSeriesSource;
    tSeriesSource = tSeriesSource - trends*b;
    %---------------------------------
    % Multitaper spectrum (need Chronux toolbox)
    %---------------------------------
    [Ss,f,Serr] = mtspectrumc(tSeriesSource,...
        struct('tapers', [7 15],'pad', 0,'fpass', [0 0.3],'err', [1,0.05],'Fs',[1/1.5]));
    %--------------------------------------------------------------------------
    % Apply filters   
    % HP
    %tSeriesSource = highpass(tSeriesSource);
    % LP
    %tSeriesSource = lowpass(tSeriesSource);
    
    
    %% Load target time series
    [c, roiInd] = intersectCols(view.coords, view.ROIs(2).coords);
    tSeriesTarget  = loadtSeries(view, 1);
    tSeriesTarget = tSeriesTarget(:,roiInd(:));
    tSeriesTarget  = raw2pc(tSeriesTarget);
    %--------------------------------------------------------------------------
    % Detrend
    %--------------------------------------------------------------------------
    trends  = rmMakeTrends(params);
    trends = single(trends);
    b = pinv(trends)*tSeriesTarget;
    tSeriesTarget = tSeriesTarget - trends*b;
    %---------------------------------
    % Multitaper spectrum (need Chronux toolbox)
    %---------------------------------
    [St,f,Serr] = mtspectrumc(tSeriesTarget,...
        struct('tapers', [7 15],'pad', 0,'fpass', [0 0.3],'err', [1,0.05],'Fs',[1/1.5]));
    %--------------------------------------------------------------------------
    % Filtering
    %--------------------------------------------------------------------------
    %HP
    %tSeriesTarget = highpass(tSeriesTarget);
    % LP
    %tSeriesTarget = lowpass(tSeriesTarget);
    
    %--------------------------------------------------------------------------
    % Plot normalized spectrum density
    %--------------------------------------------------------------------------
    figure,
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 500, 200]); %<- Set size
    plot(f,Ss/max(Ss),'Color',[255/255,204/255,0],'LineWidth', 2);
    hold on
    plot(f,St/max(St),'Color',[153/255, 204/255, 255/255],'LineWidth', 2);
    xlabel('Frequency (Hz)','FontSize',24);
    ylabel('Power (norm)','FontSize',24);
    legend('V1','V3','Location','northeast');
    legend('boxoff');
    set(gca, 'FontSize',24,'LineWidth', 2);
    set(gcf, 'color', 'w');
    set(gca, 'box', 'off');
    xlim([0 0.3]);
    filename = ['./Analysis/pwrSpectrum_'  dataTYPES(1,curDataType).name '_' roiList{1} '_' roiList{2} '.pdf'];
    export_fig(filename);
    close;
     
    %--------------------------------------------------------------------------
    % Compute CF models
    %--------------------------------------------------------------------------
    cf = computeCF(view,1,0,firstDatatype,lastDatatype,pRFmodel,VEthr);
    
    %--------------------------------------------------------------------------
    % Save time series, prf Ecc,Pol,Sigma,VE and ROI brain coordinates (3D)
    %--------------------------------------------------------------------------
    fprintf('Saving results\n');
    pRFandCFdata = struct('sourceTS',tSeriesSource,'targetTS',tSeriesTarget,...
        'sourceVE',sourceVE,'targetVE',targetVE,...
        'sourceSigma',sourceSigma,'targetSigma',targetSigma,...
        'sourceVCoords',sourceVCoords,'targetVCoords',targetVCoords,...
        'sourceCoords',sourceCoords,'targetCoords',targetCoords,...
        'sourceDistances',sourceDistances,'targetDistances',targetDistances,...
        'cf',cf);
    file = ['./Analysis/pRFandCFdata_'  dataTYPES(1,curDataType).name '_' roiList{1} '_' roiList{2} '.mat'];    
    save (file,'pRFandCFdata','-mat');
    
end

return

%---------------------------------
% High pass filter
%---------------------------------
function tSeries = highpass(tSeries)
TR = 1.5;
fc = 0.02;
fs = 1/TR;
[b, a] = butter(4,fc/(fs/2),'high');
tSeries = filtfilt(b,a,double(tSeries));
return

%---------------------------------
% Low pass filter
%---------------------------------
function tSeries = lowpass(tSeries)
TR = 1.5;
fc = 0.06;
fs = 1/TR;
[b, a] = butter(4,fc/(fs/2),'low');
tSeries = filtfilt(b,a,double(tSeries));
return

%---------------------------------
% Convert raw data to %BOLD
%---------------------------------
function data=raw2pc(data)
dc   = ones(size(data,1),1)*mean(data);
data = ((data./dc) - 1) .*100;
return;
%---------------------------------



