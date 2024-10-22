function val = sessionGet(s,param,varargin)
% Get values from a mrSESSION structure (mrVista)
%
%   val = sessionGet(s,param,varargin)
%
% This function along with sessionSet and sessionCreate should become the
% way we interact with the mrSESSION structure.  There is a long way to go
% to get these functions deeply inserted into the mrVista routines.
%
% For the new mrVista2 code, use sessGet/sessSet instead.
%
% See also: CreateNewSession (should become obsolete), sessionCreate,
% sessionSet, saveSession (should become obsolete), sessionSave (not yet
% written); 
%
% Examples:
%   s = mrSESSION;
%   v = sessionGet(s,'functionals');
%   v = sessionGet(s,'alignment');
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Dump of mrSESSION example:
%
% mrLoadRetVersion: 3.0100
%                title: ''
%             subject: 'pn'
%              examNum: []
%             inplanes: [1x1 struct]
%          functionals: [1x10 struct]
%       screenSaveSize: [512 512]
%            alignment: [4x4 double]
%          sessionCode: 'sCode'
%          description: 'sDesc'
%

global HOMEDIR

if notDefined('s'), error('mrSESSION variable required'); end
if notDefined('param'), error('Parameter field required.'); end

% make lower case and remove spaces
param = mrvParamFormat(param);

val = [];
switch param
    case 'title'
        if isfield(s, 'title'), val = s.title; end
    case 'subject'
        if isfield(s, 'subject'), val = s.subject; end
    case 'examnum'
        if isfield(s, 'examNum'), val = s.examNum; end
    case 'screensavesize'
        if isfield(s, 'screenSaveSize'), val = s.screenSaveSize; end
    case 'alignment'
        if isfield(s,'alignment'), val = s.alignment; end
    case 'sessioncode'
        if isfield(s,'sessionCode'),  val = s.sessionCode; end;
    case 'description'
        if isfield(s,'description'),val = s.description; end
        % Information about the functional scans
    case {'functionals','functionalparameters'}
        % sessionGet(s,'functionals',3);  % Third scan parameters
        % sessionGet(s,'functionals');    % All
        if isfield(s, 'functionals')
            if isempty(varargin), val = s.functionals; 
            else val = s.functionals(varargin{1});
            end
        end
    case {'pfilenames','pfilenamecellarray'}
        nScans = length(s.functionals);
        val = cell(1, nScans);
        for ii=1:nScans
            val{ii} = s.functionals(ii).PfileName;
        end
    case {'pfilelist'}
        % sessionGet(s,'pFileList',{'name1','pFile2.mag'})
        % Return indices into the functional scans corresponding to the
        % cell array of pFile names
        if isempty(varargin), error('PFile names required');  end
        allPFiles = sessionGet(s,'pFileNames');
        pFiles = varargin{1};
        val = cell(1, length(varargin));
        for ii=1:length(varargin)
            val(ii) = cellfind(allPFiles,pFiles{ii});
        end
    case {'inplanepath'}
        % pth = sessionGet(s, 'inplane path');
        % Return the path to the raw files (e.g., dicoms or nifti) for the
        % inplane anatomy (underlay for the functional data)
        
        % if the session was initialized with mrInit2, then the path to the
        % inplane should be stored in a file called mrInit2_params.mat,
        % created at initializion.
        if exist('mrInit2_params.mat', 'file') 
            p = load('mrInit2_params.mat');
            val = p.params.inplane;
        else 
            % if mrInit2_params.mat doesn't exist, try to guess the
            % location of the inplane files
           val = getIfileNames(fullfile(HOMEDIR,'Raw','Anatomy','Inplane','I'));
           if isempty(val), % perhaps the dicoms do not begin with 'I'
               val = getIfileNames(fullfile(HOMEDIR,'Raw','Anatomy','Inplane','*.dcm'));
           end           
        end
        
        if isempty(val), warning('Inplane path not found'); end %#ok<WNTAG>
        
        % if it is a char str (and not a cell array), convert to cell
        if ~isempty(val) && ~iscell(val), val = {val}; end                
        
    case {'sliceorder','sliceordering'}
        if isempty(varargin), scan = 1;
        else scan = varargin{1};
        end
        if checkfields(s.functionals(scan),'sliceOrder')
            val = s.functionals(scan).sliceOrder;
        else
            val = [];
        end
    case {'nslices','numberslices'}
        % I don't understand slquant and slices.  Older mrSESSION files
        % don't have slquant.  Newer ones appear to use it to get the
        % number of slices.  Ask for an explanation, then put it here.
        % sessionGet(mrSESSION,'nSlices')
        % sessionGet(mrSESSION,'nSlices',2)
        if isempty(varargin), scan = 1;
        else scan = varargin{1};
        end
        if checkfields(s.functionals(scan),'reconParams','slquant')
            val = s.functionals(scan).reconParams.slquant;
        else
            % I don't understand slquant and this ... we should use one but
            % not the other.  Which one?
            val = length(s.functionals(scan).slices);
        end
    case {'refslice','referenceslice','timingreferenceslice'}
        % sessionGet(mrSESSION,'refslice',2)
        if isempty(varargin), scan = 1;
        else                  scan = varargin{1};
        end
        l   = s.functionals(scan).slices;
        val = mean(l);
    case {'nshots'}
        % sessionGet(mrSESSION,'nShots',2)
        if isempty(varargin), scan = 1; 
        else                  scan = varargin{1};
        end
        if checkfields(s.functionals(scan),'reconParams','nshots')
            val = s.functionals(scan).reconParams.nshots;
        else
            val = 1;
        end
    case {'interframetiming','interframedelta','timebetweenframes','framedt'}
        % This is the proportion of a TR that separates each frame
        % acquisition. This is NOT a real number in seconds.
        % sessionGet(mrSESSION,'interframedelta',2)
        if isempty(varargin), scan = 1;
        else                  scan = varargin{1};
        end
        % When we use multiple shots, the ordering is 
        % [1a 2a 3a ... 1b 2b 3b ...] where a and b are two shots.
        % So the effective spacing between the acquisition times for slice
        % 1 and 2 is reduced by the number of shots.
        val = 1 / sessionGet(s,'nslices',scan) / sessionGet(s,'nShots',scan);
    case  {'nsamples','nframes'}
        % sessionGet(mrSESSION,'nframes',2)
        if isempty(varargin), scan = 1;
        else                  scan = varargin{1};
        end
        val = s.functionals(scan).nFrames;
    case {'tr','frameperiod'}
        % sessionGet(mrSESSION,'TR',2)
        if isempty(varargin), scan = 1;
        else                  scan = varargin{1};
        end
        val = s.functionals(scan).framePeriod;
    % Time series processing parameters for block and event analyses
    case {'eventdetrend'}
        if checkfields(s,'event','detrendFrames')
            val = s.event.detrendFrames;
        end
        
    otherwise
        error('Unknown parameter %s\n',param);
end
return;
