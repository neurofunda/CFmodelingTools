function distances = ccmVoxToVoxDistCenters(ret,rs,ROI, view, mmPerPix)
% ccmVoxToVoxDist - find distance between each pair of voxels in ROI
%     
% 2009: KVH adapted from RoiToRoiDist.

% argument checks
mrGlobals;
if ieNotDefined('view'), view = getCurView; end
if ieNotDefined('ROI') 
    if isempty(view.ROIs(1))
        error ('Must have a target ROI defined'); 
    else
        ROI = view.ROIs(1); 
    end
end 
if ieNotDefined('mmPerPix') 
    mmPerPix = readVolAnatHeader(vANATOMYPATH); 
end
if isnumeric(ROI) 
    ROI = view.ROIs(ROI); 
end
if ~strcmp(view.viewType,'Gray')
    error('Need gray-view'); 
end

% user feedback
fprintf(1,'[%s]:Computing voxel-to-voxel distances [%s]: ... ', ...
    mfilename, ROI.name); drawnow; 
tic;
% Get all nodes and edges from VOLUME view
% nodes = double(view.nodes);
% edges = double(view.edges);
r=centerCoords(1,:)
t=centerCoords(100,:)
nodes=double(r);
edges=double(t);

% Get nearest gray nodes
remappedNodes = [nodes(2,:); nodes(1,:); nodes(3,:)];
% remappedNodes(1,nodes(6,:)~=1) = 99999;
[sourceNodeIndices,dist] = nearpoints(double(ROI.coords), ...
    double(remappedNodes)); %#ok<NASGU>
[~, numSourceVoxels] = size(ROI.coords);

if length(sourceNodeIndices) ~= numSourceVoxels
    error('No gray node found for %d voxels', ...
        numSourceVoxels - length(sourceNodeIndices));
end

[lineNodeIndices,dist]=nearpoints(double(ROI.coords), ...
    double(remappedNodes)); %#ok<NASGU>
[~, numTargetVoxels] = size(ROI.coords);

if length(lineNodeIndices) ~= numTargetVoxels
    error('No gray node found for %d voxels', ...
    numTargetVoxels - length(lineNodeIndices));
end
    
% Get the distances
nodes(4,nodes(6,:)~=1) = 0;  
allDist = zeros(numel(sourceNodeIndices), numel(lineNodeIndices));


% for n=1:numel(lineNodeIndices)
%      tmp = mrManDist(nodes, edges, lineNodeIndices(n), mmPerPix, 100000, 0);
%      allDist(:,n) = tmp(sourceNodeIndices);
% end
for n=1:numel(lineNodeIndices)
     tmp = mrManDist(nodes, edges, lineNodeIndices(n), mmPerPix, 100000, 0);
     allDist(:,n) = tmp(sourceNodeIndices);
end

distances = single(allDist);


% user feedback
fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now)); drawnow;

return