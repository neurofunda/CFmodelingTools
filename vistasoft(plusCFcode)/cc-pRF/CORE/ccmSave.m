function pathStr= ccmSave(view,model,params,forceSave)
% ccmSave - save/reshape model analysis in output that mrVISTA can read
%
% 2009: KVH adapted from rmSave.

% argument checks
if notDefined('view'), error('Need view structure'); end;
if notDefined('model'), error('Need model structure'); end;
if notDefined('params'), error('Need params structure'); end;
if notDefined('forceSave'), forceSave = false; end;
    
coords = ccmGet(model{1},'roicoords');
coordsInd = ccmGet(model{1},'roiIndex');
if ~isempty(coords),
    allcoords = viewGet(view,'coords');
end;

fnames = {'x0','y0','z0','s_major','b','rss','rawrss'};
for m = 1:length(model),
    for f = 1:length(fnames),
        param = ccmGet(model{m},fnames{f});
        if numel(param) > 1 && isnumeric(param),
            if ~isempty(coords),
                switch lower(fnames{f})
                    case 'b',
                        try
                            out = zeros(1,size(allcoords,2),size(param,3));
                            out(1,coordsInd,:) = param;
                            param = out;
                        catch ME
                            fprintf(1,'[%s]: failed to get %s.\n',mfilename, fnames{f});
                            rethrow(ME);
                        end
                    otherwise
                        out = zeros(1,size(allcoords,2));
                        if numel(param)==numel(coordsInd),
                            out(coordsInd) = param;
                        else
                            out = param;
                        end
                        param = out;
                end;
            end;
            model{m} = ccmSet(model{m},fnames{f},param);
        end;
    end;
end;

params.matFileName{1} = sprintf('ccModel-%s-%s-%s.mat', ...
    view.ROIs(params.sourceROI).name, ...
    view.ROIs(params.targetROI).name, ...
    datestr(now,'yyyymmddHHMMSS'));
pathStr = fullfile(dataDir(view),params.matFileName{end});

% overwrite?
if exist(pathStr,'file') && forceSave == 0,
    [f,p] = uiputfile('*.mat','File exists, please select file?', ...
        pathStr);
    % check
    if(isequal(f,0)||isequal(p,0))
        fprintf(1,'[%s]:Model not saved.\n',mfilename);
        return;
    else
        pathStr = fullfile(p,f);
    end;
end;

params.sourceROI = view.ROIs(params.sourceROI).name;
params.targetROI = view.ROIs(params.targetROI).name;

% save
save(pathStr,'model','params');
fprintf(1,'[%s]:Saved %s.\n',mfilename,pathStr);

return;


