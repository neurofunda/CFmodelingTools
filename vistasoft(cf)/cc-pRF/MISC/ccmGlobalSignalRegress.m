function tSeries = ccmGlobalSignalRegress(view, ds, tSeries, calcPC)
% Correct for fluctuations in whole-brain signal intensity using Global
% Signal Regression (Fox et al., J.Neurophysiol., 2009).
%
% KVH 2011/04: wrote it.

gts = loadtSeries(view,ds);

if calcPC
    gts = raw2pc(gts);
end

keep = find(~isnan(sum(gts)));

% compute global mean time course
m = size(gts(:,keep),2);
g = (1/m)*gts(:,keep)*ones(m,1);

% voxel-wise regression of the time series on the global mean time course
r = pinv(g)*tSeries;

% regress the global signal out
tSeries = tSeries-(g*r);

return

function data=raw2pc(data)
    dc   = ones(size(data,1),1)*mean(data);
    data = ((data./dc) - 1) .*100;
return;