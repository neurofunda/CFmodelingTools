function [cf] = computeCF(view,nsurr,win,firstDatatype,lastDatatype,pRFmodel,threshold)
%view=initHiddenGray();
%mrGlobals;
%view = deleteAllROIs(view); view = refreshScreen(view,0);
%% Loop trough datatypes
% estimate CF's
global dataTYPES;
d=0; % datatype index for VE cell
for j=firstDatatype:lastDatatype
    %% Load datatypes
    d=d+1;
    view = selectDataType(view,j);
    curDataType = viewGet(view,'curDataType');
    %view = refreshScreen(view);
    fprintf('Loading Datatype:%s\n',dataTYPES(1,curDataType).name);
%     %% Loading ROIs
%     ROIS{1}= source; % Source
%     ROIS{2}= target; % Target
%     view = loadROI(view, ROIS,[], [], 0, 1);
%     %view = refreshScreen(view);
%     fprintf('Loading source ROI:%s\n',ROIS{1});
%     fprintf('Loading target ROI:%s\n',ROIS{2});
    %% Run ccM on original data
    surr=0;  % 0 ---> normal analysis, 1--> surrogate
    view=ccmMain_mod(view,1,2,surr,win);
    %view=ccmMain(view,1,2,surr);
    %% Loading VE of the models
    view = ccmSelect(view,2,'newest');
    view = ccmLoadDefault(view);
    % gather targetROI referred sampling size
    targetROI = viewGet(view, 'CurrentROI');
    targetCoords = view.ROIs(targetROI).coords;
    [c, roiInd] = intersectCols(view.coords, targetCoords);
    ccm = view.ccm.models{view.ccm.modelNum};
    sigma = ccm.sigma.major(roiInd);
    rss = ccm.rss(roiInd);
    rawrss = ccm.rawrss(roiInd);
    beta = ccm.beta(roiInd);
    npoints = ccm.npoints;
    ntrends = ccm.ntrends;
    %ntrends =1;
    VE=1 - (rss ./ rawrss);
    correctedVE = VE-(1-VE)*(ntrends/(npoints-ntrends-1));
    rm=load(pRFmodel);
    % check whether the source and target ROIs are loaded
    roiList=viewGet(view, 'roinames');
    sourceROI = find(strcmp(roiList,view.ccm.params.sourceROI));
    targetROI = find(strcmp(roiList,view.ccm.params.targetROI));
    % get the coordinate indices of the target ROI
    targetCoords = view.ROIs(targetROI).coords;
    [c, tRoiInd] = intersectCols(view.coords, targetCoords);
    % get all center coordinates
    centerCoords = [round(ccm.x0(tRoiInd))' ...
        round(ccm.y0(tRoiInd))' ...
        round(ccm.z0(tRoiInd))'];
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
    end
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
    
    % Following code fixes "non singleton" error. It can only handle one error only!
    % I have not found 2 in the same dataset yet but it can happen.    
    % If centerCoords = [ 0 0 0 ], set VE to zero and fix CF ecc and pol
    for i=1:length(centerCoords)
        [~,I] = intersect(centerCoords(i,:),[0 0 0],'rows');
        if I == 1
            ind(i) = 1;
        else
            ind(i) = 0;
        end
    end
    indErr = find(ind == 1);    
    %set VE to zero to bad coordinate
    VE(indErr) = 0;
    % fix CF vector length for ecc and pol
    ecc(1:indErr-1) = ecc(1:indErr-1); 
    ecc(indErr) = 0;
    ecc(indErr+1:length(VE)) = pol(indErr:length(VE)-1);
    pol(1:indErr-1) = pol(1:indErr-1); 
    pol(indErr) = 0;
    pol(indErr+1:length(VE)) = pol(indErr:length(VE)-1);
    
    %% get CF model coordinates
    CFcoords = centerCoords';     
    % Following fixes non-singleton error caused by [0 0 0] coords
    % Doesn't affect results as VE as set to 0
    CFcoords(:,indErr) = CFcoords(:,indErr-1); 
    con = zeros(length(CFcoords),1);
    cSource = zeros(length(CFcoords),3);
    indxSource = zeros(length(CFcoords),1);
    
    %% Transform coordinates into source-ROI referenced indices (provides CF indices)
    for i=1:length(CFcoords)
        [cSource(i,:), indxSource(i)] = intersectCols(coords, CFcoords(:,i));
        con(i) = indxSource(i);
    end
    
    %% Plot 3D maps
    % label = [roiList{1} '_' roiList{2}];
    %     % get old field parameters
    %     fieldNames = {'map','ph'};
    %     paramNames = {'eccentricity','polar-angle'};
    %
    %     for selfield = 1:2
    %         oldparam = viewGet(view,fieldNames{selfield});
    %         if isempty(oldparam),
    %             oldparam = cell(1,viewGet(view,'numscans'));
    %             oldparam{viewGet(view,'curscan')} = zeros(1,size(view.coords,2));
    %         end;
    %         switch lower(fieldNames{selfield})
    %             case {'map'}
    %                 oldparam{viewGet(view,'curscan')}(:) = 0;
    %                 oldparam{viewGet(view,'curscan')}(tRoiInd(keep)) = ecc;
    %                 view = viewSet(view,fieldNames{selfield},oldparam);
    %                 view = viewSet(view, 'mapName', 'eccentricity');
    %                 view = viewSet(view, 'mapUnits', char(176));
    %                 view = setClipMode(view, 'map', [0 rm.params.analysis.maxRF]);
    %             case {'ph'}
    %                 oldparam{viewGet(view,'curscan')}(:) = 0;
    %                 oldparam{viewGet(view,'curscan')}(tRoiInd(keep)) = pol;
    %                 view = viewSet(view,fieldNames{selfield},oldparam);
    %                 view = viewSet(view, 'mapName', 'polar-angle');
    %                 view = viewSet(view, 'mapUnits', 'rad');
    %             otherwise,  % do nothing
    %         end
    %
    %         % also set the colorbar title to be appropriate
    %         if checkfields(view, 'ui', 'colorbarHandle')
    %             hTitle = get(view.ui.colorbarHandle, 'Title');
    %             set(hTitle, 'String', paramNames{selfield});
    %         end
    %
    %         % refresh
    %         view  = setDisplayMode(view, fieldNames{selfield});
    %     end;
    %
    %     %% Print VE map
    %     %view = setCothresh(view, threshold);
    %     view = setDisplayMode(view,'co');
    %     updateGlobal(view);
    %     view.ui.showROIs=0;updateGlobal(view);
    %     view = meshColorOverlay(view.mesh);
    %     veMap = ['cfMap(ve)_'  dataTYPES(1,curDataType).name '_' label];
    %     fprintf('Saving eccentricity CF map\n');
    %     img = mrmGet( viewGet(view, 'Mesh'), 'screenshot' ) ./ 255; hTmp = figure('Color', 'w'); imagesc(img); axis image; axis off;   clear hTmp;
    %     export_fig(veMap);
    %     close;
    %
    %     %% Print Eccentricity map
    %     view = setCothresh(view, threshold);
    %     view = setDisplayMode(view,'map');
    %     updateGlobal(view);
    %     view.ui.showROIs=0;updateGlobal(view);
    %     view = meshColorOverlay(view.mesh);
    %     eccMap = ['cfMap(ecc)_'  dataTYPES(1,curDataType).name '_' label];
    %     fprintf('Saving eccentricity CF map\n');
    %     img = mrmGet( viewGet(view, 'Mesh'), 'screenshot' ) ./ 255; hTmp = figure('Color', 'w'); imagesc(img); axis image; axis off;  clear hTmp;
    %     export_fig(eccMap)
    %     close;
    %
    %     %% Print Polar map
    %     view = setDisplayMode(view,'ph');
    %     updateGlobal(view);
    %     view.ui.showROIs=0;
    %     view =cmapImportModeInformation(view, 'phMode', 'WedgeMapLeft_pRF.mat');
    %     updateGlobal(view);
    %     view = meshColorOverlay(view.mesh);
    %     polMap = ['cfMap(pol)_'  dataTYPES(1,curDataType).name '_' label];
    %     fprintf('Saving eccentricity CF map\n');
    %     img = mrmGet( viewGet(view, 'Mesh'), 'screenshot' ) ./ 255; hTmp = figure('Color', 'w'); imagesc(img); axis image; axis off;clear hTmp;
    %     export_fig(polMap);
    %     close;
    %
    %     %% Print CF size map
    %     view = setDisplayMode(view,'amp');
    %     updateGlobal(view);
    %     view.ui.showROIs=0;
    %     view = setClipMode(view,'amp',[0 5]);
    %     updateGlobal(view);
    %     view = meshColorOverlay(view.mesh);
    %     ampMap = ['cfMap(amp)_'  dataTYPES(1,curDataType).name '_' label];
    %     fprintf('Saving eccentricity CF map\n');
    %     img = mrmGet( viewGet(view, 'Mesh'), 'screenshot' ) ./ 255; hTmp = figure('Color', 'w'); imagesc(img); axis image; axis off;  clear hTmp;
    %     export_fig(ampMap);
    %     close;
    
    
    %% Run ccM on surrogate data
    for i=1:nsurr
        surr=1;  % 0 ---> normal analysis, 1--> surrogate
        view=ccmMain_mod(view,1,2,surr,win);
        %view=ccmMain_(view,1,2,surr,win);
        %% Loading VE of the surrogate models
        view = ccmSelect(view,2,'newest');
        view = ccmLoadDefault(view);
        % gather targetROI referred sampling size
        targetROI = viewGet(view, 'CurrentROI');
        targetCoords = view.ROIs(targetROI).coords;
        [c, roiInd] = intersectCols(view.coords, targetCoords);
        ccm = view.ccm.models{view.ccm.modelNum};
        sSigma = ccm.sigma.major(roiInd);
        sRss = ccm.rss(roiInd);
        sRawrss = ccm.rawrss(roiInd);
        sBeta = ccm.beta(roiInd);
        sVE=1 - (sRss ./ sRawrss);
        scorrectedVE = sVE-(1-sVE)*(ntrends/(npoints-ntrends-1));
        surrVE{d,i}=sVE;
        % check whether the source and target ROIs are loaded
        roiList=viewGet(view, 'roinames');
        sourceROI = find(strcmp(roiList,view.ccm.params.sourceROI));
        targetROI = find(strcmp(roiList,view.ccm.params.targetROI));
        % get the coordinate indices of the target ROI
        stargetCoords = view.ROIs(targetROI).coords;
        [c, tRoiInd] = intersectCols(view.coords, targetCoords);
        % get all center coordinates
        scenterCoords = [round(ccm.x0(tRoiInd))' ...
            round(ccm.y0(tRoiInd))' ...
            round(ccm.z0(tRoiInd))'];
        
        %         % get old field parameters
        %         fieldNames = {'map','ph'};
        %         paramNames = {'eccentricity','polar-angle'};
        %
        %         for selfield = 1:2
        %             oldparam = viewGet(view,fieldNames{selfield});
        %             if isempty(oldparam),
        %                 oldparam = cell(1,viewGet(view,'numscans'));
        %                 oldparam{viewGet(view,'curscan')} = zeros(1,size(view.coords,2));
        %             end;
        %             switch lower(fieldNames{selfield})
        %                 case {'map'}
        %                     oldparam{viewGet(view,'curscan')}(:) = 0;
        %                     oldparam{viewGet(view,'curscan')}(tRoiInd(keep)) = ecc;
        %                     view = viewSet(view,fieldNames{selfield},oldparam);
        %                     view = viewSet(view, 'mapName', 'eccentricity');
        %                     view = viewSet(view, 'mapUnits', char(176));
        %                     view = setClipMode(view, 'map', [0 rm.params.analysis.maxRF]);
        %                 case {'ph'}
        %                     oldparam{viewGet(view,'curscan')}(:) = 0;
        %                     oldparam{viewGet(view,'curscan')}(tRoiInd(keep)) = pol;
        %                     view = viewSet(view,fieldNames{selfield},oldparam);
        %                     view = viewSet(view, 'mapName', 'polar-angle');
        %                     view = viewSet(view, 'mapUnits', 'rad');
        %                 otherwise,  % do nothing
        %             end
        %
        %             % also set the colorbar title to be appropriate
        %             if checkfields(view, 'ui', 'colorbarHandle')
        %                 hTitle = get(view.ui.colorbarHandle, 'Title');
        %                 set(hTitle, 'String', paramNames{selfield});
        %             end
        %
        %             % refresh
        %             view  = setDisplayMode(view, fieldNames{selfield});
        %         end;
        %
        %         %% Print VE map
        %         %view = setCothresh(view, threshold);
        %         view = setDisplayMode(view,'co');
        %         updateGlobal(view);
        %         view.ui.showROIs=0;updateGlobal(view);
        %         view = meshColorOverlay(view.mesh);
        %         veMap = ['cfMapS(ve)_'  dataTYPES(1,curDataType).name '_' label];
        %         fprintf('Saving eccentricity CF map\n');
        %         img = mrmGet( viewGet(view, 'Mesh'), 'screenshot' ) ./ 255; hTmp = figure('Color', 'w'); imagesc(img); axis image; axis off; clear hTmp;
        %         export_fig(veMap);
        %         close;
        %
        %         %% Print Eccentricity map
        %         view = setCothresh(view, threshold);
        %         view = setDisplayMode(view,'map');
        %         updateGlobal(view);
        %         view.ui.showROIs=0;updateGlobal(view);
        %         view = meshColorOverlay(view.mesh);
        %         eccMap = ['cfMapS(ecc)_'  dataTYPES(1,curDataType).name '_' label];
        %         fprintf('Saving eccentricity CF map\n');
        %         img = mrmGet( viewGet(view, 'Mesh'), 'screenshot' ) ./ 255; hTmp = figure('Color', 'w'); imagesc(img); axis image; axis off;  clear hTmp;
        %         export_fig(eccMap)
        %         close;
        %
        %         %% Print Polar map
        %         view = setDisplayMode(view,'ph');
        %         updateGlobal(view);
        %         view.ui.showROIs=0;
        %         view =cmapImportModeInformation(view, 'phMode', 'WedgeMapLeft_pRF.mat');
        %         updateGlobal(view);
        %         view = meshColorOverlay(view.mesh);
        %         polMap = ['cfMapS(pol)_'  dataTYPES(1,curDataType).name '_' label];
        %         fprintf('Saving eccentricity CF map\n');
        %         img = mrmGet( viewGet(view, 'Mesh'), 'screenshot' ) ./ 255; hTmp = figure('Color', 'w'); imagesc(img); axis image; axis off;  clear hTmp;
        %         export_fig(polMap);
        %         close;
        %
        %         %% Print CF size map
        %         view = setDisplayMode(view,'amp');
        %         updateGlobal(view);
        %         view.ui.showROIs=0;
        %         view = setClipMode(view,'amp',[0 5]);
        %         updateGlobal(view);
        %         view = meshColorOverlay(view.mesh);
        %         ampMap = ['cfMapS(amp)_'  dataTYPES(1,curDataType).name '_' label];
        %         fprintf('Saving eccentricity CF map\n');
        %         img = mrmGet( viewGet(view, 'Mesh'), 'screenshot' ) ./ 255; hTmp = figure('Color', 'w'); imagesc(img); axis image; axis off;  clear hTmp;
        %         export_fig(ampMap);
        %         close;
        
        fprintf('surrogate # :%d\n',i);
    end
    %% Refresh ROIS
    view=deleteAllROIs(view);
    %% Save results
    sVEdist= cell2mat(surrVE);
    cf = struct('VE', VE,'correctedVE',correctedVE,'sigma',sigma,'rss',rss,'rawrss',rawrss,'beta',beta,'npoints',npoints,'ntrends',ntrends,...
        'sVE',sVEdist,'scorrectedVE',scorrectedVE,'sSigma',sSigma,'sRss',sRss,'sRawrss',sRawrss,'sBeta', sBeta,...
        'centerCoords',centerCoords,'targetCoords',targetCoords,'x0_sourceROI',x0_sourceROI,'y0_sourceROI',y0_sourceROI,'ecc',ecc,'pol',pol,...
        'scenterCoords',scenterCoords,'stargetCoords',stargetCoords,'conIndex',con);
       
%     fileID=['cfData_'  dataTYPES(1,curDataType).name '_' label '.mat'];
%     fprintf('Saving results\n');
%     save (fileID,'cf','-mat');
    
end
end




