function V=mrLoadRet3StatsToAnalyze(outFileName,scanNum,mapName)
% 
% V=mrLoadRet3StatsToAnalyze(outFileName,scanNum,mapName)
% 
% AUTHOR:  Wade
% DATE: 06.10.03
% PURPOSE: 
% Converts a statistical map in a mrLoadRet 3 Gray view into a 8 bit
% analyze file. Data are scaled from 0 to 2^8-1
% mapName can be 'co', 'amp', 'ph' or any other similar statistical map in
% If no mapName is passed, it defaults to writing out 'co'.
% If embedDimensions are passed, the data are embedded into an even larger
% dealing with BV volume anatomies that have been embedded into 256x256x256
%
% RETURNS: 
% Number of bytes written in main data block.
% 
% EXAMPLE:
% 1) V=mrLoadRet3StatsToAnalyze('temp');
% 2) V=mrLoadRet3StatsToAnalyze('temp_ph','ph');
%
% NOTES
% Originall coded to transfer mrLoadRet gray stats into SSI's EMSE / mrViewer
% $Author: wade $
% $Date: 2003/09/09 21:18:59 $
if (isempty(selectedVOLUME))
    error('You must select (click within) a volume window before proceeding');
end

if ((~exist('scanNum','var')) | (isempty(scanNum)))
% Get a coherence threshold
coThresh=getCothresh(VOLUME{selectedVOLUME});


if ((~exist('outFileName','var')) | (isempty(outFileName)))
    mapName='co';
end
volAnatSize=size(VOLUME{selectedVOLUME}.anat);
if ((~exist('VOLUME','var')) | (isempty(VOLUME{selectedVOLUME}.co)))
    error('The anatomy must be loaded');
end

% if(find(bbSize<=0))

if (~strcmp(mapName,'co'))
    error('Only co maps allowed for now');
end


% Now threshold according to the current cothresh
co(co<coThresh)=0;


blurKern=ones(3,3,3);
blurKern=blurKern./sum(blurKern(:));
blurKern=blurKern*20;

co=convn(co,blurKern,'same');



% Scale co from 0 to 2^8-1
max(co);
min(co);
% co=log10(co*10);
co(co<=0)=0;
co=co*((2^8)-1);



max(co)
min(co)

dataVolume=flipdim(dataVolume,1);
% Call SPM writevol routine to write out the data.


    s=spm_hwrite(outFileName,[ySiz xSiz zSiz],[1 1 1],1,spm_type('uint8'),0);
    V=spm_vol(outFileName);
      
    
    V.descrip=['Converted from tSeries file in session',mrSESSION.sessionCode,' : ',mrSESSION.subject,':  on ',datestr(now)];
    
    s=spm_write_vol(V,double(dataVolume));


