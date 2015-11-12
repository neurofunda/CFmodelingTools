function model = ccmSliceSet(model,tmp,slice)
% ccmSliceSet - put slice info into model struct
% we convert back to double precision here.
%
% 2009 KVH: adapted from rmSliceSet.

if ~exist('slice','var') || isempty(slice), slice = 1; end

% loop over models
for n=1:numel(model),
    % variables may have slightly different names
    ftmp   = {'x0','y0','z0','s_major','rss','rawrss'};
    % now get values from model and put in new slice values
    for fn = 1:numel(ftmp),
        val = ccmGet(model{n},ftmp{fn});
        if ~isempty(val),
            if size(val,2) == size(tmp{n}.(ftmp{fn}),2)
                val(slice,:) = double(tmp{n}.(ftmp{fn}));
            elseif slice == 1 && size(val,2) ~= size(tmp{n}.(ftmp{fn}),2)
                    val = double(tmp{n}.(ftmp{fn}));
            end
            model{n} = ccmSet(model{n},ftmp{fn},val);
        end;
    end

    % other params
    if isfield(tmp{n},'desc')
        model{n} = ccmSet(model{n},'desc',tmp{n}.desc);
    end

    % distribute beta values
    val = ccmGet(model{n},'b');
    val = val(:,:,1:size(tmp{n}.b,1));     
    for fn = 1:size(val,3),
        val(slice,:,fn) = double(tmp{n}.b(fn,:));
    end;
    model{n} = ccmSet(model{n},'b',val);
end;

return;