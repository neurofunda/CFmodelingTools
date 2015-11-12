function Ridge = ccmFit_core_ridge(data,params,dataSource)
% ccmFit_core - core of cc-pRF fit
%
% 2009: KVH wrote it.

% user feedback
 fprintf(1,'[%s]:Ridge regression: ... ', ...
     mfilename); drawnow;

tic;

% turn off warnings
warning off all


% start fitting
% for n=1:size(prediction,2)


for n=1:length(data(2,:))
    
    % minimum RSS fit
    %X = [prediction(:,n) t.trends];
    %[b,stdx,rss] = lscov(X,data);
    
    % ridge regression
    k=1;
    b = ridge(data(:,n),dataSource,k);
  

    % new rss
    rss = data(:,n)-dataSource*b;
    
    R = data(:,n)'*data(:,n) + k*diag(var(data(:,n)));
    
    % Compute RSS only for positive fits.
    nkeep = b(1,:) < 0;
    
    % Compute RSS only for positive and negative fits.
    %nkeep = b(1,:) ~= 0;
    
    rss(nkeep) = inf('single');
    
    Ridge.b{n} = b;
    Ridge.rss{n} = rss;  
    Ridge.R{n} = R;
    
    % store data with lower rss
    %minRssIndex = rss < model.rss; % this is how it was
    %     minRssIndex = rss < rss(:,n);
    
    
    % update
    %     model.x0(minRssIndex) = params.analysis.x(n);
    %     model.y0(minRssIndex) = params.analysis.y(n);
    %     model.z0(minRssIndex) = params.analysis.z(n);
    %     model.s_major(minRssIndex) = params.analysis.sigmaMajor(n);
    %     model.rss(minRssIndex) = rss(minRssIndex);
    %     %model.b([1 t.dcid+1],minRssIndex) = b(:,minRssIndex);
    
 
end;

% Correct lscov. It returns the mean rss. To maintain compatibility with
% the sum rss this function expects, we have to multiply by the divisor.
% model.rss=single(model.rss.*(size(prediction,1)-size(t.trends,2)+1));

% end time monitor
et  = toc;

% user feedback
if floor(et/3600)>0,
    fprintf(1,'Done[%d hours].\t(%s)\n', ceil(et/3600), datestr(now));
else
    fprintf(1,'Done[%d minutes].\t(%s)\n', ceil(et/60), datestr(now));
end;
drawnow;

return;


