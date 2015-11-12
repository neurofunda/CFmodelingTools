function tmp = ccmSliceGet(model,slice,id)
% ccmSliceGet - extract slice from model struct and place it
% into temporary struct, also convert to single here.
%
% 2009 KVH: adapted from rmSliceGet.

if ~exist('model','var') || isempty(model), error('Need model'); end
if ~exist('slice','var') || isempty(slice), slice = 1;           end
if ~exist('id','var') || isempty(id),       id = 1:numel(model); end

% loop over models
tmp = cell(numel(id),1);
for n=id,
    f = {'x0','y0','z0','s_major','rss','rawrss'};
    % for all models
    tmp{n}.desc = ccmGet(model{n},'desc');
    for fn = 1:numel(f),
        val = ccmGet(model{n},f{fn});
        if ~isempty(val),
            tmp{n}.(f{fn}) = single(val(slice,:));
        end;
    end;
    % put all beta values in one matrix
    val = ccmGet(model{n},'b');
    tmp{n}.b = zeros(size(val,3),size(val,2),'single');
    for fn = 1:size(val,3),
        tmp{n}.b(fn,:) = single(val(slice,:,fn));
    end;
end;

return;
