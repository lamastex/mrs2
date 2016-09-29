% script to plot one 2d function using real ranges
% plotting the 2d box and real range as a 3-d box
% 

clear variables
clear functions

%boxesFileName = 'datasets/rosenbrock_d_2_n_E4_pcf.txt';
%boxesFileName = 'datasets/sim_3_2000_pcf.txt';%'datasets/dp.txt_pcf.txt';
boxesFileName = 'datasets/dp.txt_pcf.txt';
%boxesFileName = 'datasets/sim_2_100000_pcf.txt';
%put the base name of the output files here
%this can include a path, in windows format, ie '\' path\file
outname = '';%;boxesFileName, '_Pict'];


yl=[];
% set ylim if necessary
%yl=[0 0.5];

%change the figure handle if necessary
figure;
zl=[];
% set zlim if necessary
%zl=[0 1.0];


f1 = @Function2DBoxesPlot;
    
fcol = [.8 .8 .95];

axisColour = [.3 .3 .3];

fnDiv = regexp(boxesFileName, '\.txt', 'split');
fn = fnDiv{1};
fn = regexprep(fn, '\.', 'pt');
outputfile = strcat(outname, fn, '.png');


h1 = gca;
cla(h1);

p = f1(boxesFileName, h1, fcol);
color_edge = [0.65 0.65 0.65];    % edge colour
set(p,'EdgeColor', color_edge);
%set(h1,'Zlim',[0 0.2]);
%set(h1,'Ylim',[0.5 1]);
set(h1, 'TickDir', 'out', 'XColor', axisColour, 'YColor', axisColour);
set(h1,'FontSize',14, 'FontWeight','demi');
viewangle = [21 30];
set(h1, 'View', viewangle);
if (size(zl,2) > 0) 
    set(h1,'ZLim',zl);
end;

print ('-dpng', outputfile);
