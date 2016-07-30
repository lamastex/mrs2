% txt file importing using tab delimiters
% for output from Levy example average
% shows 2 d range with range enclosure
% 

clear all

funcName = 'Levy';

boxesFileName = 'AverageHistogram.txt';
%boxesFileName = [funcName, 'Output.txt'];

dataR = dlmread(boxesFileName, '\t',0,1); % miss out the first column

base = zeros(1,size(dataR,1));

Z1 = base;
Z2 = dataR(:,2); % Collator just has straight bin heights, not counts

X1 = dataR(:,3);
X2 = dataR(:,4);

Y1 = dataR(:,5);
Y2 = dataR(:,6);

boxes = size(X1,1);

Vert=[X1(1) Y1(1) Z1(1);... % 1
    X2(1) Y1(1) Z1(1);...   % 2
    X2(1) Y2(1) Z1(1);...   % 3
    X1(1) Y2(1) Z1(1);...   % 4
    X1(1) Y1(1) Z2(1);...   % 5
    X2(1) Y1(1) Z2(1);...   % 6
    X2(1) Y2(1) Z2(1);...   % 7
    X1(1) Y2(1) Z2(1)];     % 8

FaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
Face = FaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:boxes
    new = [X1(i) Y1(i) Z1(i);... % 1
    X2(i) Y1(i) Z1(i);...   % 2
    X2(i) Y2(i) Z1(i);...   % 3
    X1(i) Y2(i) Z1(i);...   % 4
    X1(i) Y1(i) Z2(i);...   % 5
    X2(i) Y1(i) Z2(i);...   % 6
    X2(i) Y2(i) Z2(i);...   % 7
    X1(i) Y2(i) Z2(i)];     % 8

    Vert = [Vert; new];

    f = FaceBase+(8*(i-1));
    Face = [Face; f];

    t=tcolbase;
    tcol = [tcol; t];

end


boxestitle = ['Histogram for', ' ', funcName, ' average'];

figure(1);
patch('Faces',Face,'Vertices',Vert,'FaceVertexCData',tcol,...
      'FaceColor','flat')
alpha(0.3) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(boxestitle, 'FontSize', 10, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
%zlim([0 200]);




