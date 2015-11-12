function view = ccmMain(view, sourceROI, targetROI)
% ccmMain - main cc-pRF modeling program    
%
% 2009: KVH wrote it

% argument checks
if notDefined('view'), error('Need view struct'); end;

roilist = viewGet(view, 'roinames');

if isempty(roilist) 
    fprintf(1,'[%s]: Please load some ROIs first.\n', ...
            mfilename); drawnow; return
end

if ieNotDefined('sourceROI') || isempty(sourceROI)
    sourceROI = find(ccmButtonDLG('Select a source ROI',roilist));
    if isempty(sourceROI), fprintf(1,'[%s]: no source ROI selected.\n', ...
            mfilename); drawnow; return
    end
end 

if ieNotDefined('targetROI') || isempty(targetROI)
    targetROI = find(ccmButtonDLG('Select a target ROI',roilist));
    if isempty(targetROI), fprintf(1,'[%s]: no target ROI selected.\n', ...
            mfilename); drawnow; return
    end
end

view = viewSet(view,'curroi',targetROI);

% define parameters structure used in the rest of the program

params = ccmDefineParameters(view,sourceROI,targetROI);  

% user feedback:
fprintf(1,'[%s]:source ROI: [%s]; target ROI: [%s].\n',mfilename, ...
    view.ROIs(params.sourceROI).name, view.ROIs(params.targetROI).name);
fprintf(1,'[%s]:Loading data from source ROI [%s].\n',mfilename, ...
    view.ROIs(params.sourceROI).name);

% start fitting
ccmFit(view,params);
return;
