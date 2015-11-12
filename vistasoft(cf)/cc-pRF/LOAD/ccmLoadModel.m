function view = ccmLoadModel(view, varargin)
% ccmLoadModel - load data maps from cortico-cortical model file into 
% mrVista interface fields.
%
% 2009 KVH: adapted from rmLoad.

% argument checks
if ~exist('view','var') || isempty(view), view = getCurView; end;

% load file and store it in the view struct so if we load
% other parameters we should be faster.
try 
    model = view.ccm.models; 
catch %#ok<CTCH>
    model = []; 
end; 

if isempty(model),
	load(view.ccm.filename,'model');
	view.ccm.model = model;
end;

% get model names
modelNames = cell(numel(model),1);
for n=1:numel(model),
	modelNames{n} = ccmGet(model{n},'desc');
end;

% define params, but only for the ones most used
% format: 1. name for interface, 2. name for model
paramNames = {...
    'variance explained','varexplained';...
	'coherence','coherence';...
    'eccentricity','eccentricity';...
	'polar-angle','polar-angle';...
	'cc-pRF size','sigma';...
	'fit residuals (rms)','rms';...
	'x0 (mm)','x';...
	'y0 (mm)','y';...
    'z0 (mm)','z';...
	'cc-pRF betas','beta';...
    'volume','volume'};

allNamesUI = paramNames(:,1);
allNamesModel = paramNames(:,2);

% define field names
fieldNames = {'co','map','ph','amp'};

% allow for manual selection of fields:
if length(varargin) >= 3
	% we can omit the model #, but need the other two
	if isempty(varargin{1})
		sel.model = 1;
	else
		sel.model = varargin{1};
	end
	
	sel.parameter = varargin{2};
	sel.field = cellfind(fieldNames, varargin{3});
end

% get user selection if needed
if ~exist('sel','var') || isempty(sel),
    sel = rmLoadInterface(modelNames, allNamesUI, fieldNames); drawnow;
    if isempty(sel),
        fprintf('[%s]:user cancelled',mfilename);
        return;
    end
    sel.parameterName = allNamesUI{sel.parameter};
    sel.parameter = allNamesModel{sel.parameter};
else
    % input check
    if isempty(sel.model), error('No model number specified.'); end
    if isempty(sel.parameter), error('Improper parameter specification.'); 
    else sel.parameterName = sel.parameter; end
    if isempty(sel.field), error('Invalid view field specified.'); end
end

% load the variable
switch sel.parameter,
	case 'beta'
		param = ccmGet(model{sel.model},sel.parameter);
		param = param(:,:,1); % the other betas are for trends 
    otherwise
        param = ccmGet(model{sel.model},sel.parameter);
end;

if isempty(param),
	fprintf('[%s]: %s: parameter not defined.',...
		mfilename,ccmGet(model{sel.model},'desc'));
	return;
else
	% must do some sanity checks here
	param(param == Inf)  = max(param(isfinite(param(:))));
	param(param == -Inf) = min(param(isfinite(param(:))));
	param(isnan(param)) = 0;

	% if we put this in the 'co' field than we have to make sure the
	% data range from 0 to 1;
	if strcmp(fieldNames{sel.field},'co'),
        param = max(param,0);
        param = min(param,1);
	end;
    
    % if we put this in the 'ph' field than we have to make sure the
	% data range from 0 to 2*pi;
	if strcmp(fieldNames{sel.field},'ph'),
        param = max(param,0);
        param = min(param,2*pi);
	end;
end;

% get old field parameters
oldparam = viewGet(view,fieldNames{sel.field});
if isempty(oldparam),
	oldparam = cell(1,viewGet(view,'numscans'));
end;
oldparam{viewGet(view,'curscan')} = param;
view  = viewSet(view,fieldNames{sel.field},oldparam);
switch lower(fieldNames{sel.field}),
	case 'map',
		view  = viewSet(view, 'mapName', sel.parameterName);
		% also set 'ph' if you set 'co' Surface painting expects some
		% values here
		
		% set the map units for certain map types
		switch lower(sel.parameter)
			case {'eccentricity' 'sigma'}
				if ispc
					% can use special characters
					view = viewSet(view, 'mapUnits', char(176));
				else
					view = viewSet(view, 'mapUnits', char(176));
				end
			case {'polar-angle'}
				view = viewSet(view, 'mapUnits', 'rad');
			otherwise,  % do nothing
		end
		
	case 'co',
		if isempty(viewGet(view,'ph')),
			phparam = cell(1,viewGet(view,'numscans'));
			phparam{viewGet(view,'curscan')} = zeros(size(param));
			view  = viewSet(view,'ph',phparam);
		end;

		% also set the colorbar title to be appropriate
		if checkfields(view, 'ui', 'colorbarHandle')
			hTitle = get(view.ui.colorbarHandle, 'Title');
			set(hTitle, 'String', sel.parameterName);
		end
end;

% refresh
view  = setDisplayMode(view, fieldNames{sel.field});

return;
