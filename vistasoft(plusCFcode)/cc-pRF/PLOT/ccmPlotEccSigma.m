function data = ccmPlotEccSigma(vw, rm, pf, binsize, range)
% Plots cc-pRF size against pRF eccentricity
%
% KVH: adapted from rmPlotEccSigma

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

roi = vw.selectedROI;
if isequal(roi,0),
    fprintf(1,'[%s]: please load and select a ROI for plotting\n', ...
        mfilename); drawnow; return;
end

if notDefined('rm')
    rmFile = getPathStrDialog(dataDir(vw),...
        'Choose pRF model file name', 'retModel*.mat');    
    rm = load(rmFile);
end

if notDefined('binsize')
    prompt = {'Enter binsize:'};
    dlg_title = 'Please provide binsize';
    num_lines = 1;
    def = {'0.1','hsv'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    binsize = eval(answer{1});
end

if notDefined('range')
    prompt = {'Enter eccentricity range:'};
    dlg_title = 'Please provide eccentricity range';
    num_lines = 1;
    def = {['[0 ' num2str(max([rm.params.stim(:).stimSize])) ']'],'hsv'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    range = eval(answer{1});
end

if notDefined('pf'), pf = 1; end

targetROI = viewGet(vw, 'CurrentROI');

% gather targetROI referred sampling size 
targetCoords = vw.ROIs(targetROI).coords;
[c, roiInd] = intersectCols(vw.coords, targetCoords);
sig = ccm.sigma.major(roiInd);

% gather targetROI referred x0(deg) and y0(deg) 
x = rm.model{1}.x0(roiInd);
y = rm.model{1}.y0(roiInd);

%calculate eccentricity
[pol, ecc] = cart2pol(x, y);

ssr = ccm.rss(roiInd);
sst = ccm.rawrss(roiInd);
varexp = 1 - (ssr ./ sst);
varexp(~isfinite(varexp)) = 0;
varexp = max(varexp, 0);
varexp = min(varexp, 1); 

% thresholds
thresh.varexp = max(viewGet(vw,'cothresh'), rm.params.analysis.fmins.vethresh);
thresh.ecc = range + [binsize -binsize];  

% plotting parameters
xaxislim = [0 max([rm.params.stim(:).stimSize])];
MarkerSize = 8;

% find useful data given thresholds
ii = varexp > thresh.varexp & ecc > thresh.ecc(1) & ecc < thresh.ecc(2) & sig > 0.0001; 
if ~any(ii), ii = ecc > thresh.ecc(1) & ecc < thresh.ecc(2) & varexp > thresh.varexp; end

% weighted linear regression:
p = linreg(ecc(ii),sig(ii),varexp(ii));
p = flipud(p(:)); 
xfit = thresh.ecc;
yfit = polyval(p,xfit);

% output struct
data.ecc = ecc(ii);
data.sig = sig(ii);
data.ve  = varexp(ii);
data.xfit = xfit(:);
data.yfit = yfit(:);
data.x    = (thresh.ecc(1):binsize:thresh.ecc(2))';
data.y    = nan(size(data.x));
data.ysterr = nan(size(data.x));

% compute the 95% bootstrap confidence interval if possible 
if exist('bootstrp','file')
    B = bootstrp(1000,@(x) localfit(x,ecc(ii),sig(ii),varexp(ii)),(1:numel(ecc(ii))));
    B = B';
    pct1 = 100*0.05/2;
    pct2 = 100-pct1;
    b_lower = prctile(B',pct1);
    b_upper = prctile(B',pct2);
    if pf
        fprintf(1,'[%s]: intercept = [%0.2f:%0.2f]; slope = [%0.2f:%0.2f]\n', ...
            mfilename, b_lower(1), b_upper(1), b_lower(2), b_upper(2));
    end
    keep1 = B(1,:)>b_lower(1) &  B(1,:)<b_upper(1);
    keep2 = B(2,:)>b_lower(2) &  B(2,:)<b_upper(2);
    keep = keep1 & keep2;
    data.b_xfit = linspace(min(xfit),max(xfit),100)';
    fits = [ones(100,1) data.b_xfit]*B(:,keep);
    data.b_upper = max(fits,[],2);
    data.b_lower = min(fits,[],2);  
end 

% plot averaged data
for b=thresh.ecc(1):binsize:thresh.ecc(2),
    bii = ecc >  b-binsize./2 & ...
          ecc <= b+binsize./2 & ii;
    if any(bii),
        s = wstat(sig(bii),varexp(bii));
        % ii2 = find(data.x==b);
        ii2 = abs(data.x - b) < 0.000001;
        data.y(ii2) = s.mean;
        data.ysterr(ii2) = s.sterr;       
    end;
end;

% plot if requested
if pf
    figure('Color', 'w'); hold on;
    
    if ~exist('bootstrp','file')
        errorbar(data.x,data.y,data.ysterr,'ko',...
            'MarkerFaceColor','k',...
            'MarkerSize',MarkerSize);
    else
        p1 = plot(data.b_xfit,data.b_upper,'--','LineWidth',1.5);
        p2 = plot(data.b_xfit,data.b_lower,'--','LineWidth',1.5);
        set(p1,'Color',[.6 .6 .6]);
        set(p2,'Color',[.6 .6 .6]);               
    end
    
    plot(data.x,data.y,'ko',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',MarkerSize);
    
    plot(xfit,yfit','k','LineWidth',2);
    ylabel('Connective Field Radius (mm)');xlabel('Eccentricity (deg)');
    % axis([xaxislim(1) xaxislim(2) 0 max(data.y)+1]);
    axis([xaxislim(1) xaxislim(2) 0 15]);
    title([vw.ccm.params.sourceROI ' >> ' vw.ROIs(targetROI).name]);
end

return

function  B=localfit(ii,x,y,ve)
    B = linreg(x(ii),y(ii),ve(ii));
    B(:);
return
