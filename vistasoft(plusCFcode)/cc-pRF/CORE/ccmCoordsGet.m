function val = ccmCoordsGet(viewType, model, param, coords)
% ccmCoordsGet - wrapper for ccmGet to limit values to certain coordinates
%
% KVH: adapted from rmCoordsGet.

if ~exist('viewType','var') || isempty(viewType),   error('Need view type'); end;
if ~exist('model','var')    || isempty(model),      error('Need model');     end;
if ~exist('param','var')    || isempty(param),      error('Need param');     end;
if ~exist('coords','var'),                          error('Need coords');    end;
if isempty(coords), val = []; return; end
% allow the model to be a cell array -- just take the first entry
if iscell(model), model = model{1};  end

tmp = ccmGet(model, param);

switch lower(viewType),
    case 'inplane'
        val = zeros(size(coords, 2), 1);
        for n = 1:length(val),
            
            val(n) = tmp(coords(1,n), coords(2,n), coords(3,n));
        end;

    case 'gray'
        % allow 3xN gray coords to be specified, as well as gray node
        % indices:
        if size(coords, 1)==3 
            % 3xN coords specification: remap into indices
            allCoords = viewGet(getSelectedGray, 'coords');
            [c, coords] = intersectCols(allCoords, coords);
        end
        
        if numel(size(tmp)) == 2,
            val = tmp(coords);
		else
			if length(coords)==1
				% the squeeze command will mis-orient fields like beta, by
				% permuting across multiple dimensions. We want it to
				% voxels x predictors.
				val = permute(tmp(1,coords,:), [2 3 1]);
			else
	            val = squeeze(tmp(1,coords,:));
			end
		end
		
	otherwise, error('Invalid view type.')
end;

return
