function ccmMenu = ccmAnalysisMenu(analysismenu, vw)
% Submenu for cortico-cortical pRF modeling
%
% 2012 KVH: wrote it.

ccmMenu = uimenu(analysismenu,'Label','Cortico-Cortical Model', ...
    'Separator','off');

cb=[vw.name,'=ccmMain(',vw.name,',[],[]);', ...
    vw.name,'=refreshScreen(',vw.name,');'];
uimenu(ccmMenu,'Label','Run (cc-pRF)', ...
    'Separator','off', 'CallBack', cb);  

return