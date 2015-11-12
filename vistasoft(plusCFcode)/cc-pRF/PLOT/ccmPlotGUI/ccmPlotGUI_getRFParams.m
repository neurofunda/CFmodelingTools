function rfParams = ccmPlotGUI_getRFParams(model, viewType, coords)
% For the cortico-cortical model GUI: get cc-pRF params for a single voxel.
%
% 2011 KVH: adapted from rmPlotGUI_getRFparams

rfParams = zeros(1,4);
        
rfParams(1) = ccmCoordsGet(viewType, model, 'x0', coords);
rfParams(2) = ccmCoordsGet(viewType, model, 'y0', coords);
rfParams(3) = ccmCoordsGet(viewType, model, 'z0', coords);
rfParams(4) = ccmCoordsGet(viewType, model, 'sigmamajor',coords);
        
end
