function ccmPlotGUIsw(vw, roi, voxel, update)
% rmPlotGUI - create a GUI for visualizing the results of a
% cortico-cortical analysis
%
% 2011 KVH: adapted from rmPlotGUI.
% 2014 NG: modified it

if notDefined('vw') || isempty('vw'), vw = getCurView; end
if notDefined('update'), update = 0; end

% switch: allow the GUI controls to call this function requesting an update
if update, 
      ccmPlotGUIsw_update(vw); 
      %ccmPlotGUIswVFM_update(vw); 
    return;          
end;

try
	ccm = vw.ccm.models;
catch %#ok<CTCH>
    ccm = [];
end;

if isempty(ccm)
    fprintf(1,'[%s]: please load a cortico-cortical model first\n', ...
        mfilename); drawnow; return;
end

roiList=viewGet(vw, 'roinames');
sourceROI = find(strcmp(roiList,vw.ccm.params.sourceROI), 1);
if isempty(sourceROI)
    fprintf(1,'[%s]: please load %s first\n',mfilename,...
        vw.ccm.params.sourceROI); 
    drawnow; return;
end

if ~exist('roi','var') || isempty(roi),
    roi = vw.selectedROI;
    if isequal(roi,0), 
        fprintf(1,'[%s]: please load and select a ROI for plotting\n', ...
            mfilename); drawnow; return;
    end
end;

% disambiguate ROI specification: name, ROI index, coords, struct...
% this will return an ROI struct
roi = tc_roiStruct(vw, roi);

% grab information from the model and ROI into a compact description 
% (the struct variable M). M will be stored as the UserData property of the
% GUI figure.
M = ccmPlotGUI_getModel_mod(vw, roi);

% open the GUI
M = ccmPlotGUIsw_openFig(M,voxel);

% run an initial refresh
ccmPlotGUIsw_update(vw, M);
%ccmPlotGUIswVFM_update(vw, M);

return;
