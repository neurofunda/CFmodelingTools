function ccmMenu = ccmLoadMenu(fileMenu, view)
% Submenu for cortico-cortical pRF modeling
%
% 2012 KVH: wrote it.

ccmMenu = uimenu(fileMenu,'Label','Cortico-Cortical Model','Separator','off');

% Load cc-pRF model with default parameters
callBackstr=[view.name ' = ccmSelect(',view.name,', 2); ' ...
		     view.name ' = ccmLoadDefault(' view.name '); '];
uimenu(ccmMenu, 'Label', 'Select and load cortico-cortical model', 'Separator', 'off',...
		'CallBack', callBackstr);

% Load cc-pRF model derived eccentricity in interface
callBackstr=[view.name ' = ccmLoadDerived(' view.name '); '...
		     view.name ' = refreshScreen(' view.name '); '];
uimenu(ccmMenu,'Label','Load derived maps','Separator','off',...
    'CallBack',callBackstr);

return