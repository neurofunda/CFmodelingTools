function ccmPlotGUI(view, roi)
% rmPlotGUI - create a GUI for visualizing the results of a
% cortico-cortical analysis
%
% 2011 KVH: adapted from rmPlotGUI.

if notDefined('view'), view = getCurView; end

% switch: allow the GUI controls to call this function requesting an update
if isequal(view, 'update'), 
    ccmPlotGUI_update; 
    return;          
end;

try
	ccm = view.ccm.models;
catch %#ok<CTCH>
    ccm = [];
end;

if isempty(ccm)
    fprintf(1,'[%s]: please load a cortico-cortical model first\n', ...
        mfilename); drawnow; return;
end

roiList=viewGet(view, 'roinames');
sourceROI = find(strcmp(roiList,view.ccm.params.sourceROI), 1);
if isempty(sourceROI)
    fprintf(1,'[%s]: please load %s first\n',mfilename,...
        view.ccm.params.sourceROI); 
    drawnow; return;
end

if ~exist('roi','var') || isempty(roi),
    roi = view.selectedROI;
    if isequal(roi,0), 
        fprintf(1,'[%s]: please load and select a ROI for plotting\n', ...
            mfilename); drawnow; return;
    end
end;

% disambiguate ROI specification: name, ROI index, coords, struct...
% this will return an ROI struct
roi = tc_roiStruct(view, roi);

% grab information from the model and ROI into a compact description 
% (the struct variable M). M will be stored as the UserData property of the
% GUI figure.
M = ccmPlotGUI_getModel(view, roi);

% open the GUI
M = ccmPlotGUI_openFig(M);

% run an initial refresh
ccmPlotGUI_update(M);

return;
