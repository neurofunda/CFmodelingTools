function view = ccmLoadDerived(view, varargin)
% ccmLoadDerived - derives retinotopic maps using known map in source ROI 
% and a cortico-cortical model into mrVista interface fields
%
% 2011 KVH: wrote it

% argument checks
if ~exist('view','var') || isempty(view), view = getCurView; end;

try
	ccm = view.ccm.models{view.ccm.modelNum};
catch %#ok<CTCH>
    ccm = [];
end;

if isempty(ccm)
    fprintf(1,'[%s]: please load a cortico-cortical model first\n', ...
        mfilename); drawnow; return;
end

% Load the retinotopic model file
rmFile = getPathStrDialog(dataDir(view),...
        'Choose pRF model file name', 'retModel*.mat');    
rm = load(rmFile);

% check whether the source and target ROIs are loaded
roiList=viewGet(view, 'roinames');
sourceROI = find(strcmp(roiList,view.ccm.params.sourceROI));
targetROI = find(strcmp(roiList,view.ccm.params.targetROI));
if isempty(sourceROI) && isempty(targetROI)
    disp(['Please first load ' view.ccm.params.sourceROI ...
        ' and ' view.ccm.params.targetROI]); return    
elseif isempty(sourceROI)
    disp(['Please first load ' view.ccm.params.sourceROI]); return
elseif isempty(targetROI)
    disp(['Please first load ' view.ccm.params.targetROI]); return
end
viewSet(view,'curroi',targetROI);

% get the coordinate indices of the target ROI
targetCoords = view.ROIs(targetROI).coords;
[c, tRoiInd] = intersectCols(view.coords, targetCoords);

% get all center coordinates
centerCoords = [round(ccm.x0(tRoiInd))' ...
                round(ccm.y0(tRoiInd))' ...
                round(ccm.z0(tRoiInd))'];
            
% now find corresponding coordIndices
fprintf(1,'[%s]: Deriving the retinotopic map in %s from %s ... ', ...
    mfilename, view.ROIs(targetROI).name, view.ROIs(sourceROI).name);
drawnow; h = waitbar(0,'Please wait...'); tic
coords = view.ROIs(sourceROI).coords;
[coords roiInd] = intersectCols(view.coords, coords);
x0_sourceROI = []; y0_sourceROI = []; keep = [];
for i = 1:length(centerCoords)
    for j = 1:length(coords)
        if centerCoords(i,:)' == coords(:,j) %#ok<BDSCA>
            x0_sourceROI = [x0_sourceROI rm.model{1}.x0(roiInd(j))]; %#ok<AGROW>
            y0_sourceROI = [y0_sourceROI rm.model{1}.y0(roiInd(j))]; %#ok<AGROW>
            keep = [keep i]; %#ok<AGROW>
            break
        end
    end
    waitbar(i / length(centerCoords));
end
et = toc; close(h); fprintf(1,'done [%0.2f sec].\n', et);

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


% get old field parameters
fieldNames = {'map','ph'};
paramNames = {'eccentricity','polar-angle'};

for selfield = 1:2
    oldparam = viewGet(view,fieldNames{selfield});
    if isempty(oldparam),
        oldparam = cell(1,viewGet(view,'numscans'));
        oldparam{viewGet(view,'curscan')} = zeros(1,size(view.coords,2));
    end;            
    switch lower(fieldNames{selfield})
        case {'map'}
            oldparam{viewGet(view,'curscan')}(:) = 0;
            oldparam{viewGet(view,'curscan')}(tRoiInd(keep)) = ecc;
            view = viewSet(view,fieldNames{selfield},oldparam);
            view = viewSet(view, 'mapName', 'eccentricity');
            view = viewSet(view, 'mapUnits', char(176));
            view = setClipMode(view, 'map', [0 rm.params.analysis.maxRF]);
        case {'ph'}
            oldparam{viewGet(view,'curscan')}(:) = 0;
            oldparam{viewGet(view,'curscan')}(tRoiInd(keep)) = pol;
            view = viewSet(view,fieldNames{selfield},oldparam);
            view = viewSet(view, 'mapName', 'polar-angle');
            view = viewSet(view, 'mapUnits', 'rad');
        otherwise,  % do nothing
    end
    
    % also set the colorbar title to be appropriate
    if checkfields(view, 'ui', 'colorbarHandle')
        hTitle = get(view.ui.colorbarHandle, 'Title');
        set(hTitle, 'String', paramNames{selfield});
    end
    
    % refresh
    view  = setDisplayMode(view, fieldNames{selfield});
end;
view.ui.mapMode = setColormap(view.ui.mapMode, 'hsvTbCmap');

return;
