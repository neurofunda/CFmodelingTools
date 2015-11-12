function [reply ok] = ccmButtonDLG(headerStr,optionStr)
% Button Dialog Box
%
%  reply = buttondlg(headerStr,optionStr,[defaultRes])
%
% Method for creating a window that allows the user to toggle and select
% options in the cell array string 'optionStr'.
%
% INPUTS:
%  
%  headerStr: string for the title of the dialog.
%
%	optionStr: cell array of strings, one for each button option.
%
%	defaultRes: logical array of default responses, 1 for selected, 0 for
%	unselected. [Default: all zeros, don't select anything]
%
% OUTPUTS:
%
% reply: a boolean vector with length(optionStr) that indicates
%        selected options. reply is empty if the user
%        chooses to cancel.
%
% ok: flag indicating whether the user canceled or not. 
%
% EXAMPLE:
%  reply = buttondlg('pick it',{'this','that','the other'})

ok = 0;
OptionsPerColumn=30; % max number of options in one column

if nargin<2
    disp('Error: "bottondlg" requires two inputs');
    return
end

if isunix,  fontSize = 10;
else        fontSize = 9;
end

if iscell(optionStr), optionStr = char(optionStr); end

nOptions = size(optionStr,1);

ncols=ceil(nOptions/OptionsPerColumn);
nOptionsPerColumn=ceil(nOptions/ncols);

% scale factors for x and y axis coordinates
xs = 1.8;  
ys = 1.8;

% default sizes
butWidth=10;
botMargin = 0.2;
height = nOptionsPerColumn+2+botMargin*2+.5;

% If we don't have a minimum colwidth (5 in this case), short strings don't
% show up.
colwidth=max(size(optionStr,2),5);%
width = max(ncols*colwidth,length(headerStr))+2;
width = max(width,2*butWidth+2);

% open the figure
h = figure('MenuBar', 'none',...
    'Units', 'char',...
    'Resize', 'on',...
    'NumberTitle', 'off',...
    'Position', [20, 10, width*xs,height*ys]); %     
% center the figure -- we needed to use char to get
% a reasonable size, but want to move the corners to 
% a centered position. This is a quick way to do it:
set(h, 'Units', 'normalized');
normPos = get(h, 'Position'); % pos in normalized units
normPos(1:2) = [.5 .5] - normPos(3:4)./2;
set(h, 'Position', normPos);


% Display title text
x = 1;
y = nOptionsPerColumn+1+botMargin;

uicontrol('Style','text',...
    'Units','char',...
    'String',headerStr,...
    'Position',[(x+.25)*xs,(y+.3)*ys,(width-2.5)*xs,ys*.9],...
    'HorizontalAlignment','center',...
    'FontSize',fontSize);

% Display the radio buttons
y = y+botMargin/2;
y0=y;
c=0;

h_buttongroup = uibuttongroup('visible','off');
for optionNum=1:nOptions
    y = y-1;
    h_button(optionNum) = uicontrol('Style','radio',...
              'Units','char',...
              'String',optionStr(optionNum,:),...
              'Position',[x*xs+c*colwidth*xs,y*ys,colwidth*xs+(c+1)*colwidth*xs,ys],...
              'HorizontalAlignment','left',...
              'parent', h_buttongroup, ...
              'FontSize',fontSize); %#ok<AGROW>
    if optionNum >= (c+1)*nOptionsPerColumn
        c=c+1;
        y=y0;
    end
end

set(h_buttongroup,'SelectedObject',[]);  % No selection
set(h_buttongroup,'Visible','on');

% Cancel button
x=1;
y=botMargin;
uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Units','char',...
    'Position',[x*xs,y*ys,butWidth*xs/2,ys],...
    'CallBack','uiresume',...
    'FontSize',fontSize,...
    'UserData','Cancel');

% OK button
x = width-butWidth-1;
uicontrol('Style','pushbutton',...
    'String','OK',...
    'Units','char',...
    'Position',[x*xs/2+2,y*ys,butWidth*xs/2,ys],...
    'CallBack','uiresume',...
    'FontSize',fontSize,...
    'UserData','OK');

% let the user select some radio buttons and
% wait for a 'uiresume' callback from OK/Cancel
uiwait

% determine which button was hit.
response = get(gco,'UserData');

% gather the status of the radio buttons if 'OK' was 
% selected.  Otherwise return empty matrix.
if strcmp(response,'OK')
    for optionNum=1:nOptions
        reply(optionNum)=get(h_button(optionNum),'Value'); %#ok<AGROW>
    end
	ok = 1;
else
    reply = [];
end

close(h)

return;

  