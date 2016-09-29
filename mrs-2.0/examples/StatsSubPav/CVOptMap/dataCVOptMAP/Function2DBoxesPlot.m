function [ patchhandle ] = Function2DBoxesPlot(inputfile, axeshandle, varargin )
%Read data from inputfile and plot onto axeshandle.

%Function is set up for 2-D data with real range.
%Plots 2-D box and real range as a 3-d box

%Input:
%inputfile:     full text filename of file where intput is found
%               eg 'myinput.txt'
%axeshandle:    handle of axes on which to put the plot
%
%Up to 4 optional input:
%               colour for top and bottom faces of the boxes plotted
%               colour for edges of the faces of the boxes plotted
%               face alpha (between 0 for clear, 1 opaque, inclusive)
%               edge alpha (between 0 for clear, 1 opaque, inclusive)
%
%Output:
%The handle of the patch graphic.

dataR = dlmread(char(inputfile), '\t', 0, 1); % from row 0, col 1

OneZ2 = dataR(:,2);
OneZ1 = zeros(size(OneZ2,1));

OneX1 = dataR(:,3);
OneX2 = dataR(:,4);

OneY1 = dataR(:,5);
OneY2 = dataR(:,6);

Oneboxes = size(OneX1,1);

OneVert=[OneX1(1) OneY1(1) OneZ1(1);... % 1
    OneX2(1) OneY1(1) OneZ1(1);...   % 2
    OneX1(1) OneY2(1) OneZ1(1);...   % 3 4
    OneX2(1) OneY2(1) OneZ1(1);...   % 4 3
    OneX1(1) OneY1(1) OneZ2(1);...   % 5
    OneX2(1) OneY1(1) OneZ2(1);...   % 6
    OneX1(1) OneY2(1) OneZ2(1);...   % 7 8
    OneX2(1) OneY2(1) OneZ2(1)];     % 8 7

OneFaceBase = [1 2 4 3;...     % bottom
    5 6 8 7;...             % top
    1 2 6 5;...             % front
    4 3 7 8;...             % back
    2 4 8 6;...             % rhs    
    3 1 5 7];               % lhs
   


OneFace = OneFaceBase;

fcol = [.7 0.6 1];
if  size(varargin,2) > 0
    fcol = varargin{1};
end

tcolbase = [fcol;...
    fcol;...
    1 1 1;...
    1 1 1;...
    1 1 1;...
    1 1 1];
tcol = tcolbase;

for i=2:Oneboxes
    new = [OneX1(i) OneY1(i) OneZ1(i);... % 1
    OneX2(i) OneY1(i) OneZ1(i);...   % 2
    OneX1(i) OneY2(i) OneZ1(i);...   % 3
    OneX2(i) OneY2(i) OneZ1(i);...   % 4
    OneX1(i) OneY1(i) OneZ2(i);...   % 5
    OneX2(i) OneY1(i) OneZ2(i);...   % 6
    OneX1(i) OneY2(i) OneZ2(i);...   % 7
    OneX2(i) OneY2(i) OneZ2(i)];     % 8

    OneVert = [OneVert; new];

    f = OneFaceBase+(8*(i-1));
    OneFace = [OneFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end

ecol = [.25 .25 .25];    % edge colour
if  size(varargin,2) > 1
    ecol = varargin{2};
end

falpha = 0.6;    % edge colour
if  size(varargin,2) > 2
    falpha = varargin{3};
end

ealpha = 1;    % edge colour
if  size(varargin,2) > 3
    ealpha = varargin{4};
end

axes(axeshandle);

patchhandle = patch('Faces',OneFace,'Vertices',OneVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','FaceAlpha', falpha,'EdgeColor', ecol,'EdgeAlpha',ealpha);
xlim([min(OneX1) max(OneX2)]);
ylim([min(OneY1) max(OneY2)]);
zlim([0 max(OneZ2)]);


end

