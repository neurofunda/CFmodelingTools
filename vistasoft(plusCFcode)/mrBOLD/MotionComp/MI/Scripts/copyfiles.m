% This script copies the datasets into a working directory and executes
%         fileList = {'at041009','at041016','ctr040710','ctr040717',...
fileList = {'rd20050427'};
% if isunix
%     networkPathTarget = '/biac2/wandell';
%     networkPathSource = '/snarp';
% else
%     networkPathTarget = '\\White\biac2-wandell\';
%     networkPathSource = '\\snarp'
% end
%   temp dir on teal
for i = 1:length(fileList)
register
%temp dir on biac2
% if isunix
copyList = dir(tempDir);
for i = 3:length(copyList)