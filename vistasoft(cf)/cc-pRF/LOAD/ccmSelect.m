function view = ccmSelect(view, loadModel, ccmFile)
% ccmSelect - select cortico-cortical model file and put in view struct
%
% 2009 KVH: adapted from rmSelect

% argument checks
if ~exist('view','var') || isempty(view), view = getCurView; end;
if ~exist('loadModel','var') || isempty(loadModel), loadModel = true; end;

% choose filename:
if ~exist('ccmFile','var') || isempty(ccmFile),
    ccmFile = getPathStrDialog(dataDir(view),...
        'Choose cortico-cortical model file name', ...
        'ccModel*.mat');
    drawnow;
elseif iscell(ccmFile)
    ccmFile = ccmFile{1};
end

% if user just wants the newest file, check for it:
if ischar(ccmFile) && ismember(lower(ccmFile), {'newest' 'mostrecent'})
	pattern = fullfile( dataDir(view), 'ccModel-*.mat' );
	w = dir(pattern);
	if isempty(w)
		error('Most Recent File selected; no ccModel-* files found.')
	end
	[dates order] = sortrows( datevec({w.date}) ); % oldest -> newest
	ccmFile = fullfile( dataDir(view), w(order(end)).name );
end

% if the load model flag is 1, but the file's already selected, just load
% it and return:
if loadModel && checkfields(view,'ccm','filename')
    if ~exist(ccmFile, 'file')
        error('File not found'); 
    end;
    load(ccmFile,'model','params');
    view.ccm.filename = ccmFile;
    view.ccm.models   = model;
    view.ccm.params   = params;
    view.ccm.modelNum = 1;
    return;
end;

if ~exist(ccmFile,'file') 
    if exist([ccmFile '.mat'], 'file')
        ccmFile = [ccmFile '.mat'];
    elseif check4File(fullfile(dataDir(view), ccmFile))
        ccmFile = fullfile(dataDir(view), ccmFile);
    else
        fprintf('[%s]:No file: %s',mfilename,ccmFile);
        return;
    end
end
    
% store ccmFile filename:
view.ccm.filename = ccmFile;

if loadModel==0
    % clear previous models but don't load them untill we need them:
    view.ccm.models = [];
else
    % go ahead and load
    load(ccmFile,'model','params');
    view.ccm.models = model;
    view.ccm.params = params;
    view.ccm.modelNum = 1;
end;
    
return;

