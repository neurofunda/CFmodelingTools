function params = ccmDefineParameters(view, sourceROI, targetROI)
% ccmDefineParameters - define parameters for cc-pRF model fit
%
% 2009 KVH: wrote it.

% argument checks
if ~exist('view','var') || isempty(view), error('Need view struct'); end
if ieNotDefined('sourceROI'), error('Need source ROI index'); end
if ieNotDefined('targetROI'), error('Need target ROI index'); end

params.sourceROI = sourceROI;
params.targetROI = targetROI;

% All layer-1 voxels in the sourceROI are cc-pRF centers
params.analysis.centers = 1:length(view.ROIs(params.sourceROI).coords);

% How many cc-pRF sizes to test? 
params.analysis.sigmas = linspace(0.0001,25,50);

% Which CC-pRF model do we want to use? So far only a circular symmetric
% Gaussian is implemented.
params.analysis.pRFmodel = {'one gaussian'};

% Convert data to percent BOLD signal?
params.analysis.calcPC = true;

% apply low-pass filtering (fc = cutoff frequency; 0 = no filtering)
params.analysis.fc = 0.1; 

% remove global signals?
params.analysis.gsr = 0;

% Having the 'nSlices' field be definable as a stimulus
% argument creates bugs -- so, we make a separate, analysis parameter
% here. This may not be a complete solution, but seems like the most
% appropriate fix; this function should be run within rmMain for each view:
params.analysis.nSlices = viewGet(view, 'numSlices');

% Keep track of some other meta-information, such
% as the scanning session and data type of the data used. This will be
% useful for finding data, in those functions where the view structure is
% not passed as an argument.
try
    mrGlobals;
    params.analysis.session = mrSESSION.sessionCode;
    params.analysis.dataType = getDataTypeName(view);
    params.analysis.scans = (1:viewGet(view, 'numScans'));
    params.analysis.viewType = view.viewType;
catch %#ok<CTCH>
    fprintf('[%s]: Couldn''t record meta-data about scan params.\n', mfilename);
end

% Gather the necessary scan parameters
params.stim = dataTYPES(view.curDataType).scanParams;

% Detrending is done using dicrete cosine transforms:
for n = 1:numel(params.stim)
    params.stim(n).nDCT = 3;
end

% There is never a repetition of the time-course in the source ROI
for n = 1:numel(params.stim)
    params.stim(n).nUniqueRep = 1;
end

% set output file name base
if ~exist('matFileName','var') || isempty(matFileName),
        params.matFileName{1} = ['ccModel-',datestr(now,'yyyymmdd-HHMMSS')];
else
    if ~iscell(matFileName),
        params.matFileName{1} = matFileName;
    else
        params.matFileName = matFileName;
    end
end;

% We could allow temporal smoothing. This can be achieved by
% decimation, allowing both SNR improvement due to low-pass filtering and
% speed up due to less sample points.
params.analysis.decimate = 0;

% tSeries predictions are based on the activity in the source ROI
params.analysis.allstimimages = ccmLoadData(view,params,params.sourceROI);

return


