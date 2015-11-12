function ccmPlotLinkage(vw, savefigures)
% Plots the pRF centers of the target ROI against the pRF centers of the
% source ROI. Also reports some stats.
%
% 2010 KVH: wrote it 

MarkerSize = 6;
nbins = 50;

% check view struct
if notDefined('vw'), vw = getCurView; end
if notDefined('savefigures'), savefigures = 0; end

% cortico-cortical model must be loaded
try
	ccm = vw.ccm.models{vw.ccm.modelNum};
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

targetROI = viewGet(vw, 'CurrentROI');

% Load the retinotopic model file
rmFile = getPathStrDialog(dataDir(vw),...
        'Choose pRF model file name', 'retModel*.mat');    
rm = load(rmFile);

targetCoords = vw.ROIs(targetROI).coords;
[c, roiInd] = intersectCols(vw.coords, targetCoords);

% compute variance explained
ssr = ccm.rss(roiInd);
sst = ccm.rawrss(roiInd);
varexp = 1 - (ssr ./ sst);
varexp(~isfinite(varexp)) = 0;
varexp = max(varexp, 0);
varexp = min(varexp, 1); 

% first get all center coordinates
centerCoords = [round(ccm.x0(roiInd))' ...
                round(ccm.y0(roiInd))' ...
                round(ccm.z0(roiInd))'];

% gather targetROI referred x0(deg) and y0(deg) 
x0_targetROI = rm.model{1}.x0(roiInd);
y0_targetROI = rm.model{1}.y0(roiInd);
                        
% now find corresponding coordIndices
fprintf(1,'[%s]: Linking pRF centers in %s with those in %s ... ', ...
    mfilename, vw.ROIs(targetROI).name, vw.ROIs(sourceROI).name);
drawnow; tic
coords = vw.ROIs(sourceROI).coords;
[coords roiInd] = intersectCols(vw.coords, coords);
x0_sourceROI = []; 
y0_sourceROI = []; keep = [];
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
et = toc; fprintf(1,'done [%0.2f sec].\n', et);

x0_targetROI = x0_targetROI(keep); 
y0_targetROI = y0_targetROI(keep);
varexp = varexp(keep);

%calculate eccentricity
[sPol, sEcc] = cart2pol(x0_sourceROI, y0_sourceROI);
[tPol, tEcc] = cart2pol(x0_targetROI, y0_targetROI);

%calculate polar angle
sPol = cart2pol(x0_sourceROI, y0_sourceROI);
sPol = mod(sPol, 2*pi);
tPol = cart2pol(x0_targetROI, y0_targetROI);
tPol = mod(tPol, 2*pi);


%% Plot the average eccentricity data

% bin size (eccentricity range) of the data
binsize = max([rm.params.stim(:).stimSize])/nbins;

thresh.varexp = viewGet(vw,'cothresh');
thresh.sEcc   = [0 max([rm.params.stim(:).stimSize])]; 
thresh.tEcc   = [0 max([rm.params.stim(:).stimSize])]; 

ii = varexp > thresh.varexp & ...
     sEcc > thresh.sEcc(1) & sEcc < thresh.sEcc(2) & ...
     tEcc > thresh.tEcc(1) & tEcc < thresh.tEcc(2);
 
% report how much of the variance the unity line explains  

x = sEcc(ii)'; 
y = tEcc(ii)';

[r,p] = corr(x,y);
fprintf(1,'[%s]: Eccentricity:\t correlation = %0.2f; p = %0.2f\n', ...
        mfilename, r, p);

% output struct
data.x      = (thresh.sEcc(1):binsize:thresh.sEcc(2))';
data.y      = nan(size(data.x));
data.ysterr = nan(size(data.x));

% plot averaged data
for b=thresh.sEcc(1):binsize:thresh.sEcc(2),
    bii = sEcc >  b-binsize./2 & ...
          sEcc <= b+binsize./2 & ii;
    if any(bii),
        % weighted mean of sigma
        s = wstat(tEcc(bii),varexp(bii));
        % store
        ii2 = find(data.x==b);
        data.y(ii2) = s.mean;
        data.ysterr(ii2) = s.sterr;
    end;
end;

hFig1 = figure('Color', 'w'); hold on;
errorbar(data.x,data.y,data.ysterr,'ko', 'MarkerFaceColor','k', 'MarkerSize', MarkerSize);
h1 = refline(1); set(h1, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
axis square
xlabel(['Eccentricity (deg) source ROI: ' vw.ROIs(sourceROI).name]);
ylabel(['Eccentricity (deg) target ROI: ' vw.ROIs(targetROI).name]);
axis([0 thresh.sEcc(2) 0 thresh.sEcc(2)]);

if savefigures, 
    saveas(hFig1,['Linkage_Eccentricity_' ...
        viewGet(vw,'datatype') '_' ...
        vw.ROIs(sourceROI).name '_' ...
        vw.ROIs(targetROI).name '.png']);
end

%% Plot the average polar angle data

binsize = (2*pi)/nbins;

thresh.sPol = [0 2*pi]; 
thresh.tPol = [0 2*pi]; 

ii = varexp > thresh.varexp & ...
     sPol > thresh.sPol(1) & sPol < thresh.sPol(2) & ...
     tPol > thresh.tPol(1) & tPol < thresh.tPol(2) & ...
     sEcc > thresh.sEcc(1) & sEcc < thresh.sEcc(2) & ...
     tEcc > thresh.tEcc(1) & tEcc < thresh.tEcc(2);

% report how much of the variance the refline explains

x = sPol(ii)'; 
y = tPol(ii)';

[rN,pN] = corr(x,y); [rC,pC] = circ_corrcc(x,y);
fprintf(1,'[%s]: Polar angle:\t correlation = %0.2f; p = %0.2f\n', ...
        mfilename, max(rN,rC), min(pN,pC));

% output struct
data.x      = (thresh.sPol(1):binsize:thresh.sPol(2))';
data.y      = nan(size(data.x));
data.ysterr = nan(size(data.x));

% plot averaged data
for b=thresh.sPol(1):binsize:thresh.sPol(2),
    bii = sPol >  b-binsize./2 & ...
          sPol <= b+binsize./2 & ii;
    if any(bii),
        % weighted mean of tPol
        s = circ_stats(tPol(bii),varexp(bii)');
        % store
        ii2 = find(data.x==b);
        data.y(ii2) = mod(s.mean,2*pi);
        data.ysterr(ii2) = s.sterr;
    end;
end;

hFig2 = figure('Color', 'w'); hold on;
errorbar(data.x,data.y,data.ysterr,'ko', 'MarkerFaceColor','k', 'MarkerSize', MarkerSize);
h1 = refline(1); set(h1, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
axis square
xlabel(['Polar angle (rad) source ROI: ' vw.ROIs(sourceROI).name]);
ylabel(['Polar angle (rad) target ROI: ' vw.ROIs(targetROI).name]);
axis([0 2*pi 0 2*pi]); 

if savefigures, 
    saveas(hFig2,['Linkage_PolarAngle_' ...
        viewGet(vw,'datatype') '_' ...
        vw.ROIs(sourceROI).name '_' ...
        vw.ROIs(targetROI).name '.png']);
end

return