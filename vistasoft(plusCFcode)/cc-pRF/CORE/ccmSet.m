function model = ccmSet(model,param,val)
% ccmSet - set data from various cortico-cortical models
%
% 2009 KVH: adapted from rmSet.

if nargin == 0,
	% initiate model, bit circular but useful
	model.description = 'init';
	model = ccmSet(model,'x0');
	model = ccmSet(model,'y0');
    model = ccmSet(model,'z0');
	return;
else
	% set a specific parameter
	if ieNotDefined('model'), error('model not defined');    end;
	if ieNotDefined('param'), error('param not defined');    end;
	if ieNotDefined('val'),   val = [];                      end;
end;

switch lower(param),
	case {'desc','description'}
		model.description = val;
	case {'x0'}
		model.x0 = val;
	case {'y0'}
		model.y0 = val;
    case {'z0'}            
        model.z0 = val;
	case {'s','sigma'}
		model.sigma.major = val;
		model.sigma.minor = val;
		model.sigma.theta = zeros(size(val));
	case {'sigmamajor','sigma major','s_major'}
		model.sigma.major = val;
	case {'sigmaminor','sigma minor','s_minor'}
		model.sigma.minor = val;
	case {'sigmatheta','sigma theta','s_theta'}
		model.sigma.theta = val;
    case {'b', 'beta'}
        model.beta = val;
    case {'nt','ntrends','number of trends'},
        model.ntrends = val;
	case {'rss','residual sum of squares'}
		model.rss = val;
	case {'rawrss','raw residual sum of squares'}
		model.rawrss = val;
	case {'coords','roicoords'}
		model.roi.coords = val;
    case {'indices','roiindices','roiind','roiindex'}
        model.roi.coordsIndex = val;
    case {'roiname','roi name'}
        model.roi.name = val;
	case {'n','npoints','number of data points'}
		model.npoints = val;
	otherwise,
		error('Unknown parameter: %s',param);
end;
