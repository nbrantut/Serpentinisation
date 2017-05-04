function exportfig(figname, varargin)
%EXPORTFIG sets up a figure for export in publication format.
%
%usage:
%
% EXPORTFIG(figname)
%
%Without options, export the current figure and generates a pdf file 
%based on figname.
%The default size is 8.4cm width by 6.5cm height, and is the standard single
%column JGR figure (in the good old format!).
%
% EXPORTFIG(figname, 'xSize', x, 'ySize', y)
%
%Using optional argument pairs 'xSize' and/or 'ySize' allows to change the
%figure size. Values must be in cm.
%
% EXPORTFIG(figname, 'hfig', h)
%
%Using optional argument pair 'hfig', h allows to use a given figure handle
%to save.
%
%input:
%   figname:    a string of character with the name of the figure.
%   'xSize', x: sets the figure width to x (in cm).
%   'ySize', y: sets the figure height to y (in cm).
%   'hfig', h:  uses the figure handle h to export.


p = inputParser;

p.addRequired('figname', @ischar);
p.addParamValue('xSize', 8.4, @(x) isnumeric(x) && x>0);
p.addParamValue('ySize', 6.5, @(x) isnumeric(x) && x>0);
p.addParamValue('font','Times',@ischar);
p.addParamValue('fontsize',9,@isnumeric);
p.addParamValue('hfig', get(0,'CurrentFigure'), @ishandle);
p.parse(figname, varargin{:}); 

xSize = p.Results.xSize;
ySize = p.Results.ySize;
font = p.Results.font;
fontsize = p.Results.fontsize;
hfig = p.Results.hfig;

% plot

%find objects which have text:
htext = findall(hfig,'-property','FontName');

for k=1:length(htext)
    set(htext(k), 'FontName', font, ...
        'FontSize',fontsize);
end

set(hfig,'PaperUnits','centimeters');

%Additional coordinates to center the figure on A4-paper
xLeft = (21.0-xSize)/2;
yTop = (29.7-ySize)/2;

set(hfig,'PaperPosition',[xLeft yTop xSize ySize])

% save
filename = strcat(figname,'.eps');

if strfind(version, '2016')
    %add .eps to make filename
    filenametmp = strcat(figname,'_tmp.eps');
    %print
    print(hfig,'-depsc',filenametmp);
    %fix hyphen
    cmd = ['sed -e ''s/\/hyphen/\/endash/g'' -e ''s/[1[[:space:]]3][[:space:]]0[[:space:]]setdash/[2 2] 0 setdash/g'' ' filenametmp ' > ' filename];
    system(cmd);
    cmd = ['rm -f ' filenametmp];
    system(cmd);
else
    %print
    print(hfig,'-depsc',filename);
    fix_dottedline(filename);
end
%close(hfig);

if ismac || isunix
    cmd = ['ps2pdf -dEPSCrop -dAutoRotatePages=/None '  filename];
    system(cmd);
    cmd = ['rm -f ' filename];
    system(cmd);
else
    eps2pdf([figname '_tmp.eps'],'C:\Program Files\gs\gs9.21\bin\gswin64c.exe')
    copyfile([figname '_tmp.pdf'],[figname '.pdf'])
    delete([figname '_tmp.pdf'])
    delete([figname '_tmp.eps'])
    delete([figname '.eps'])
end
    
    

