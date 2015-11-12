function val = ccmGet(model,param,varargin)
% ccmGet - retrieve data from various cortico-cortical models
%
% 2009 KVH: adapted from rmGet.

if ~exist('model','var') || isempty(model), error('No model defined.');    end;
if ~exist('param','var') || isempty(param), error('No parameter defined'); end;

% default output
val = [];

try
    switch lower(param),
        case {'x0'}
            val = model.x0;
        case {'y0'}
            val = model.y0;
        case {'z0'}
            val = model.z0;
        case {'desc','description'}
            val = model.description;
        case {'s','sigma'}
            val = (model.sigma.major + model.sigma.minor)./2;
        case {'sigmamajor','sigma major','s_major'}
            val = model.sigma.major;
        case {'sigmaminor','sigma minor','s_minor'}
            val = model.sigma.minor;
        case {'sigmatheta','sigma theta','s_theta'}
            val = model.sigma.theta;
        case {'b','beta'}
            val = model.beta;
        case {'nt','ntrends','number of trends'},
            %val = model.beta.ntrends;
            val = model.ntrends;
        case {'varexp','varexplained','variance','variance explained', 'varianceexplained'}
            warning('off','MATLAB:divideByZero');
            val = 1 - (model.rss ./ model.rawrss);
            val(~isfinite(val)) = 0;
            val = max(val, 0);
            val = min(val, 1);
        case {'coh','coherence'}
            warning('off','MATLAB:divideByZero');
            val = sqrt(ccmGet(model,'variance explained'));
        case {'rss','residual sum of squares'}
            val = model.rss;
        case {'rawrss','raw sum of squares'}
            val = model.rawrss;
        case {'rms'}
            val = sqrt(model.rss./model.npoints);
        case {'coords','roicoords'}
            val = model.roi.coords;
        case {'indices','roiindices','roiind','roiindex'}
            val = model.roi.coordsIndex;
        case {'roiname','roi name'},
            val = model.roi.name;
        case {'n','npoints','number of data points'}
            val = model.npoints;
        otherwise,
            error('[%s]:Unknown parameter: %s.',mfilename,param);
    end;
catch  %#ok<CTCH>
    val = [];
end

return;
