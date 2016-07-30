% txt file importing using tab delimiters
% for output from BivGaussian example two histograms
% shows 2 d range with height
%

clear all

funcName = 'BivGaussian';



figure(3);

% --------------- histogram split on k -----------------------
Oneboxestitle = 'Split on max count = 100';

boxesFileName = 'BivGaussianFirst.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out first column

Onebase = zeros(1,size(dataR,1));

OneCounts = dataR(:,2);
OneCountsTotal = sum(OneCounts);

OneZ1 = Onebase;
OneZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
OneZ2 = OneZ2/OneCountsTotal; % height as relative_count_in_box/volume

OneX1 = dataR(:,3);
OneX2 = dataR(:,4);

OneY1 = dataR(:,5);
OneY2 = dataR(:,6);

Oneboxes = size(OneX1,1);

OneVert=[OneX1(1) OneY1(1) OneZ1(1);... % 1
    OneX2(1) OneY1(1) OneZ1(1);...   % 2
    OneX2(1) OneY2(1) OneZ1(1);...   % 3
    OneX1(1) OneY2(1) OneZ1(1);...   % 4
    OneX1(1) OneY1(1) OneZ2(1);...   % 5
    OneX2(1) OneY1(1) OneZ2(1);...   % 6
    OneX2(1) OneY2(1) OneZ2(1);...   % 7
    OneX1(1) OneY2(1) OneZ2(1)];     % 8

OneFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs
    4 1 5 8];               % lhs
OneFace = OneFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Oneboxes
    new = [OneX1(i) OneY1(i) OneZ1(i);... % 1
    OneX2(i) OneY1(i) OneZ1(i);...   % 2
    OneX2(i) OneY2(i) OneZ1(i);...   % 3
    OneX1(i) OneY2(i) OneZ1(i);...   % 4
    OneX1(i) OneY1(i) OneZ2(i);...   % 5
    OneX2(i) OneY1(i) OneZ2(i);...   % 6
    OneX2(i) OneY2(i) OneZ2(i);...   % 7
    OneX1(i) OneY2(i) OneZ2(i)];     % 8

    OneVert = [OneVert; new];

    f = OneFaceBase+(8*(i-1));
    OneFace = [OneFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(1,2,1);
patch('Faces',OneFace,'Vertices',OneVert,'FaceVertexCData',tcol,...
      'FaceColor','flat')
alpha(0.3) % transparency
view(39.0, 10.0)
set(gca,'ZGrid','on')
title(Oneboxestitle, 'FontSize', 10, 'FontName', 'Ariel', 'FontWeight', 'Bold')
zlim([0 0.2]);
zlabel('proportion/vol');

% --------------- histogram with priorty queue split -----------------------
Twoboxestitle = 'Priority queue split until 50 leaves';

boxesFileName = 'BivGaussianSecond.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1);

Twobase = zeros(1,size(dataR,1));


TwoCounts = dataR(:,2);
TwoCountsTotal = sum(TwoCounts);

TwoZ1 = Twobase;
TwoZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
TwoZ2 = TwoZ2/TwoCountsTotal; % height as relative_count_in_box/volume

TwoX1 = dataR(:,3);
TwoX2 = dataR(:,4);

TwoY1 = dataR(:,5);
TwoY2 = dataR(:,6);

Twoboxes = size(TwoX1,1);

TwoVert=[TwoX1(1) TwoY1(1) TwoZ1(1);... % 1
    TwoX2(1) TwoY1(1) TwoZ1(1);...   % 2
    TwoX2(1) TwoY2(1) TwoZ1(1);...   % 3
    TwoX1(1) TwoY2(1) TwoZ1(1);...   % 4
    TwoX1(1) TwoY1(1) TwoZ2(1);...   % 5
    TwoX2(1) TwoY1(1) TwoZ2(1);...   % 6
    TwoX2(1) TwoY2(1) TwoZ2(1);...   % 7
    TwoX1(1) TwoY2(1) TwoZ2(1)];     % 8

TwoFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs
    4 1 5 8];               % lhs
TwoFace = TwoFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Twoboxes
    new = [TwoX1(i) TwoY1(i) TwoZ1(i);... % 1
    TwoX2(i) TwoY1(i) TwoZ1(i);...   % 2
    TwoX2(i) TwoY2(i) TwoZ1(i);...   % 3
    TwoX1(i) TwoY2(i) TwoZ1(i);...   % 4
    TwoX1(i) TwoY1(i) TwoZ2(i);...   % 5
    TwoX2(i) TwoY1(i) TwoZ2(i);...   % 6
    TwoX2(i) TwoY2(i) TwoZ2(i);...   % 7
    TwoX1(i) TwoY2(i) TwoZ2(i)];     % 8

    TwoVert = [TwoVert; new];

    f = TwoFaceBase+(8*(i-1));
    TwoFace = [TwoFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(1,2,2);
patch('Faces',TwoFace,'Vertices',TwoVert,'FaceVertexCData',tcol,...
      'FaceColor','flat')
alpha(0.3) % transparency
view(39.0, 10.0)
set(gca,'ZGrid','on')
title(Twoboxestitle, 'FontSize', 10, 'FontName', 'Ariel', 'FontWeight', 'Bold')
zlim([0 0.2]);
zlabel('proportion/vol');

