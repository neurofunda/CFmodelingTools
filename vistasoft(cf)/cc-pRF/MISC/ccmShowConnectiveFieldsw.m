function vw = ccmShowConnectiveFieldsw(vw,rfParams,M,distances,target,refSource)
% ccmLoadDefault: load the default maps from a cortico-cortical model.
%
% 2012 KVH: wrote it
% 2014 NG: modified it

if notDefined('vw'), vw = getCurView; end
vw = getCurView;  %% mod for show on mesh plot to work

map = cell(1,viewGet(vw,'numscans'));
map{viewGet(vw,'curscan')} = zeros(1,size(vw.coords,2));
map{viewGet(vw,'curscan')}(:) = 0;


roiList=viewGet(vw, 'roinames');
M.sourceROI = find(strcmp(roiList,vw.ccm.params.sourceROI));
M.targetROI = find(strcmp(roiList,vw.ccm.params.targetROI)); %% nuevo
M.sCoords   = vw.ROIs(M.sourceROI).coords;
M.params    = vw.ccm.params;
%M.distances = ccmVoxToVoxDist(vw.ROIs(M.sourceROI), vw, vw.mmPerVox);
M.distances = distances;

RF = zeros(numel(M.params.roi.coordsIndex), 1);

for i = 1:size(M.sCoords,2)
    if rfParams(1) == M.sCoords(1,i) && ...
            rfParams(2) == M.sCoords(2,i) && ...
            rfParams(3) == M.sCoords(3,i)
        X = M.distances(:,i);
        RF = single(exp(-1.*((X.^2)./(2.*rfParams(4).^2))));
        break;
    end
end

[~, roiInd] = intersectCols(vw.coords, M.sCoords);


%% ridge regression
% if get(M.ui.showRidge, 'Value') == 1 && any(rfParams)
%     map{viewGet(vw,'curscan')}(roiInd) = cell2mat(vw.ridge.b(M.prevVoxel))*15;
% else
    map{viewGet(vw,'curscan')}(roiInd) = RF ; %b(1,:) < 0
% end    
    vw = viewSet(vw, 'map', map);
    
    map{viewGet(vw,'curscan')}(:) = -1;
    map{viewGet(vw,'curscan')}(roiInd) = 1;
    vw = viewSet(vw, 'co', map);
    
    vw = setClipMode(vw, 'map', [0 1]); updateGlobal(vw);
    
    vw.ui.mapMode = setColormap(vw.ui.mapMode, 'hotCmap');
    
    %% mod to plot target roi also
    % load target ROI CF
    vw = makePointROI(vw,target,1); updateGlobal(vw);
    vw.ROIs(1,3).color = [1 0 0];
    vw.ui.showROIs = -2;
  
 
    updateGlobal(vw);
    vw = setDisplayMode(vw, 'map');
    vw = refreshScreen(vw);
    vw = meshColorOverlay(vw);   
    
    % delete target ROI
    vw = deleteROI(vw,3); vw=refreshScreen(vw,0); updateGlobal(vw);
   
    
    return
