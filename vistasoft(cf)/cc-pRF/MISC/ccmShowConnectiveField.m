function vw = ccmShowConnectiveField(vw,rfParams,distances)
% ccmLoadDefault: load the default maps from a cortico-cortical model.
%
% 2012 KVH: wrote it
if notDefined('vw'), vw = getCurView; end
vw = getCurView;  %% mod  for show on mesh plot to work

map = cell(1,viewGet(vw,'numscans'));
map{viewGet(vw,'curscan')} = zeros(1,size(vw.coords,2));
map{viewGet(vw,'curscan')}(:) = 0;

roiList=viewGet(vw, 'roinames');
M.sourceROI = find(strcmp(roiList,vw.ccm.params.sourceROI));
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
map{viewGet(vw,'curscan')}(roiInd) = RF;
vw = viewSet(vw, 'map', map);
 
map{viewGet(vw,'curscan')}(:) = -1;
map{viewGet(vw,'curscan')}(roiInd) = 1;
vw = viewSet(vw, 'co', map);

vw = setClipMode(vw, 'map', [0 1]);

vw.ui.mapMode = setColormap(vw.ui.mapMode, 'jetCmap');

updateGlobal(vw);
vw = setDisplayMode(vw, 'map');
vw = refreshScreen(vw);
vw = meshColorOverlay(vw);

return
