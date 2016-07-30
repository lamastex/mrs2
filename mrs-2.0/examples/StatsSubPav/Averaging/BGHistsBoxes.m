% txt file importing using tab delimiters
% for output from Averaging Bivariate Gaussian example for separate histograms
% shows 2 d range with range enclosure
% 
% 
clear all

funcName = 'BivGaussian';
%boxestitle = ['Histograms for 10', ' ', funcName, ' samples'];

%obtain density first

Mus = zeros(1,2);
CovMat = eye(2);
[SGx, SGy] = meshgrid(-5:.1:5,-5:.2:5);
SGXY = [SGx(:),SGy(:)];
SGZ = EvalGaussPDF(SGXY,Mus,CovMat);
SGz = reshape(SGZ,size(SGx));


figure(2);

% --------------- histogram number One -----------------------
Oneboxestitle = 'k = 26';

boxesFileName = 'BivGaussian1.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Onebase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

OneZ1 = Onebase;
OneZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
OneZ2 = OneZ2/CountsTotal; % height as relative_count_in_box/volume

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

color_edge = [.25 .25 .25];    % edge colour

subplot(5,2,1);
patch('Faces',OneFace,'Vertices',OneVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Oneboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.21]);
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Two -----------------------
Twoboxestitle = 'k = 52';

boxesFileName = 'BivGaussian2.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Twobase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

TwoZ1 = Twobase;
TwoZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
TwoZ2 = TwoZ2/CountsTotal; % height as relative_count_in_box/volume

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


subplot(5,2,2);
patch('Faces',TwoFace,'Vertices',TwoVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Twoboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.2]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Three -----------------------
Threeboxestitle = 'k = 78';

boxesFileName = 'BivGaussian3.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Threebase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

ThreeZ1 = Threebase;
ThreeZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
ThreeZ2 = ThreeZ2/CountsTotal; % height as relative_count_in_box/volume

ThreeX1 = dataR(:,3);
ThreeX2 = dataR(:,4);

ThreeY1 = dataR(:,5);
ThreeY2 = dataR(:,6);

Threeboxes = size(ThreeX1,1);

ThreeVert=[ThreeX1(1) ThreeY1(1) ThreeZ1(1);... % 1
    ThreeX2(1) ThreeY1(1) ThreeZ1(1);...   % 2
    ThreeX2(1) ThreeY2(1) ThreeZ1(1);...   % 3
    ThreeX1(1) ThreeY2(1) ThreeZ1(1);...   % 4
    ThreeX1(1) ThreeY1(1) ThreeZ2(1);...   % 5
    ThreeX2(1) ThreeY1(1) ThreeZ2(1);...   % 6
    ThreeX2(1) ThreeY2(1) ThreeZ2(1);...   % 7
    ThreeX1(1) ThreeY2(1) ThreeZ2(1)];     % 8

ThreeFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
ThreeFace = ThreeFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Threeboxes
    new = [ThreeX1(i) ThreeY1(i) ThreeZ1(i);... % 1
    ThreeX2(i) ThreeY1(i) ThreeZ1(i);...   % 2
    ThreeX2(i) ThreeY2(i) ThreeZ1(i);...   % 3
    ThreeX1(i) ThreeY2(i) ThreeZ1(i);...   % 4
    ThreeX1(i) ThreeY1(i) ThreeZ2(i);...   % 5
    ThreeX2(i) ThreeY1(i) ThreeZ2(i);...   % 6
    ThreeX2(i) ThreeY2(i) ThreeZ2(i);...   % 7
    ThreeX1(i) ThreeY2(i) ThreeZ2(i)];     % 8

    ThreeVert = [ThreeVert; new];

    f = ThreeFaceBase+(8*(i-1));
    ThreeFace = [ThreeFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(5,2,3);
patch('Faces',ThreeFace,'Vertices',ThreeVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Threeboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.15]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Four -----------------------
Fourboxestitle = 'k = 104';

boxesFileName = 'BivGaussian4.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Counts = dataR(:,2);
CountsTotal = sum(Counts);

Fourbase = zeros(1,size(dataR,1));

%Z1 = dataR(:,3)./dataR(:,2);  
FourZ1 = Fourbase;
FourZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
FourZ2 = FourZ2/CountsTotal; % height as relative_count_in_box/volume


FourX1 = dataR(:,3);
FourX2 = dataR(:,4);

FourY1 = dataR(:,5);
FourY2 = dataR(:,6);

Fourboxes = size(FourX1,1);

FourVert=[FourX1(1) FourY1(1) FourZ1(1);... % 1
    FourX2(1) FourY1(1) FourZ1(1);...   % 2
    FourX2(1) FourY2(1) FourZ1(1);...   % 3
    FourX1(1) FourY2(1) FourZ1(1);...   % 4
    FourX1(1) FourY1(1) FourZ2(1);...   % 5
    FourX2(1) FourY1(1) FourZ2(1);...   % 6
    FourX2(1) FourY2(1) FourZ2(1);...   % 7
    FourX1(1) FourY2(1) FourZ2(1)];     % 8

FourFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
FourFace = FourFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Fourboxes
    new = [FourX1(i) FourY1(i) FourZ1(i);... % 1
    FourX2(i) FourY1(i) FourZ1(i);...   % 2
    FourX2(i) FourY2(i) FourZ1(i);...   % 3
    FourX1(i) FourY2(i) FourZ1(i);...   % 4
    FourX1(i) FourY1(i) FourZ2(i);...   % 5
    FourX2(i) FourY1(i) FourZ2(i);...   % 6
    FourX2(i) FourY2(i) FourZ2(i);...   % 7
    FourX1(i) FourY2(i) FourZ2(i)];     % 8

    FourVert = [FourVert; new];

    f = FourFaceBase+(8*(i-1));
    FourFace = [FourFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(5,2,4);
patch('Faces',FourFace,'Vertices',FourVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Fourboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.15]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Five -----------------------
Fiveboxestitle = 'k = 130';

boxesFileName = 'BivGaussian5.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Fivebase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

%Z1 = dataR(:,3)./dataR(:,2);  
FiveZ1 = Fivebase;
FiveZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
FiveZ2 = FiveZ2/CountsTotal; % height as relative_count_in_box/volume


FiveX1 = dataR(:,3);
FiveX2 = dataR(:,4);

FiveY1 = dataR(:,5);
FiveY2 = dataR(:,6);

Fiveboxes = size(FiveX1,1);

FiveVert=[FiveX1(1) FiveY1(1) FiveZ1(1);... % 1
    FiveX2(1) FiveY1(1) FiveZ1(1);...   % 2
    FiveX2(1) FiveY2(1) FiveZ1(1);...   % 3
    FiveX1(1) FiveY2(1) FiveZ1(1);...   % 4
    FiveX1(1) FiveY1(1) FiveZ2(1);...   % 5
    FiveX2(1) FiveY1(1) FiveZ2(1);...   % 6
    FiveX2(1) FiveY2(1) FiveZ2(1);...   % 7
    FiveX1(1) FiveY2(1) FiveZ2(1)];     % 8

FiveFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
FiveFace = FiveFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Fiveboxes
    new = [FiveX1(i) FiveY1(i) FiveZ1(i);... % 1
    FiveX2(i) FiveY1(i) FiveZ1(i);...   % 2
    FiveX2(i) FiveY2(i) FiveZ1(i);...   % 3
    FiveX1(i) FiveY2(i) FiveZ1(i);...   % 4
    FiveX1(i) FiveY1(i) FiveZ2(i);...   % 5
    FiveX2(i) FiveY1(i) FiveZ2(i);...   % 6
    FiveX2(i) FiveY2(i) FiveZ2(i);...   % 7
    FiveX1(i) FiveY2(i) FiveZ2(i)];     % 8

    FiveVert = [FiveVert; new];

    f = FiveFaceBase+(8*(i-1));
    FiveFace = [FiveFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(5,2,5);
patch('Faces',FiveFace,'Vertices',FiveVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Fiveboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.125]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Six -----------------------
Sixboxestitle = 'k = 156';

boxesFileName = 'BivGaussian6.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Sixbase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

%Z1 = dataR(:,3)./dataR(:,2);  
SixZ1 = Sixbase;
SixZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
SixZ2 = SixZ2/CountsTotal; % height as relative_count_in_box/volume

SixX1 = dataR(:,3);
SixX2 = dataR(:,4);

SixY1 = dataR(:,5);
SixY2 = dataR(:,6);

Sixboxes = size(SixX1,1);

SixVert=[SixX1(1) SixY1(1) SixZ1(1);... % 1
    SixX2(1) SixY1(1) SixZ1(1);...   % 2
    SixX2(1) SixY2(1) SixZ1(1);...   % 3
    SixX1(1) SixY2(1) SixZ1(1);...   % 4
    SixX1(1) SixY1(1) SixZ2(1);...   % 5
    SixX2(1) SixY1(1) SixZ2(1);...   % 6
    SixX2(1) SixY2(1) SixZ2(1);...   % 7
    SixX1(1) SixY2(1) SixZ2(1)];     % 8

SixFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
SixFace = SixFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Sixboxes
    new = [SixX1(i) SixY1(i) SixZ1(i);... % 1
    SixX2(i) SixY1(i) SixZ1(i);...   % 2
    SixX2(i) SixY2(i) SixZ1(i);...   % 3
    SixX1(i) SixY2(i) SixZ1(i);...   % 4
    SixX1(i) SixY1(i) SixZ2(i);...   % 5
    SixX2(i) SixY1(i) SixZ2(i);...   % 6
    SixX2(i) SixY2(i) SixZ2(i);...   % 7
    SixX1(i) SixY2(i) SixZ2(i)];     % 8

    SixVert = [SixVert; new];

    f = SixFaceBase+(8*(i-1));
    SixFace = [SixFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(5,2,6);
patch('Faces',SixFace,'Vertices',SixVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Sixboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.125]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Seven -----------------------
Sevenboxestitle = 'k = 182';

boxesFileName = 'BivGaussian7.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Sevenbase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

%Z1 = dataR(:,3)./dataR(:,2);  
SevenZ1 = Sevenbase;
SevenZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
SevenZ2 = SevenZ2/CountsTotal; % height as relative_count_in_box/volume


SevenX1 = dataR(:,3);
SevenX2 = dataR(:,4);

SevenY1 = dataR(:,5);
SevenY2 = dataR(:,6);

Sevenboxes = size(SevenX1,1);

SevenVert=[SevenX1(1) SevenY1(1) SevenZ1(1);... % 1
    SevenX2(1) SevenY1(1) SevenZ1(1);...   % 2
    SevenX2(1) SevenY2(1) SevenZ1(1);...   % 3
    SevenX1(1) SevenY2(1) SevenZ1(1);...   % 4
    SevenX1(1) SevenY1(1) SevenZ2(1);...   % 5
    SevenX2(1) SevenY1(1) SevenZ2(1);...   % 6
    SevenX2(1) SevenY2(1) SevenZ2(1);...   % 7
    SevenX1(1) SevenY2(1) SevenZ2(1)];     % 8

SevenFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
SevenFace = SevenFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Sevenboxes
    new = [SevenX1(i) SevenY1(i) SevenZ1(i);... % 1
    SevenX2(i) SevenY1(i) SevenZ1(i);...   % 2
    SevenX2(i) SevenY2(i) SevenZ1(i);...   % 3
    SevenX1(i) SevenY2(i) SevenZ1(i);...   % 4
    SevenX1(i) SevenY1(i) SevenZ2(i);...   % 5
    SevenX2(i) SevenY1(i) SevenZ2(i);...   % 6
    SevenX2(i) SevenY2(i) SevenZ2(i);...   % 7
    SevenX1(i) SevenY2(i) SevenZ2(i)];     % 8

    SevenVert = [SevenVert; new];

    f = SevenFaceBase+(8*(i-1));
    SevenFace = [SevenFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(5,2,7);
patch('Faces',SevenFace,'Vertices',SevenVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Sevenboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.125]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Eight -----------------------
Eightboxestitle = 'k = 208';

boxesFileName = 'BivGaussian8.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Eightbase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

%Z1 = dataR(:,3)./dataR(:,2);  
EightZ1 = Eightbase;
EightZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
EightZ2 = EightZ2/CountsTotal; % height as relative_count_in_box/volume

EightX1 = dataR(:,3);
EightX2 = dataR(:,4);

EightY1 = dataR(:,5);
EightY2 = dataR(:,6);

Eightboxes = size(EightX1,1);

EightVert=[EightX1(1) EightY1(1) EightZ1(1);... % 1
    EightX2(1) EightY1(1) EightZ1(1);...   % 2
    EightX2(1) EightY2(1) EightZ1(1);...   % 3
    EightX1(1) EightY2(1) EightZ1(1);...   % 4
    EightX1(1) EightY1(1) EightZ2(1);...   % 5
    EightX2(1) EightY1(1) EightZ2(1);...   % 6
    EightX2(1) EightY2(1) EightZ2(1);...   % 7
    EightX1(1) EightY2(1) EightZ2(1)];     % 8

EightFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
EightFace = EightFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Eightboxes
    new = [EightX1(i) EightY1(i) EightZ1(i);... % 1
    EightX2(i) EightY1(i) EightZ1(i);...   % 2
    EightX2(i) EightY2(i) EightZ1(i);...   % 3
    EightX1(i) EightY2(i) EightZ1(i);...   % 4
    EightX1(i) EightY1(i) EightZ2(i);...   % 5
    EightX2(i) EightY1(i) EightZ2(i);...   % 6
    EightX2(i) EightY2(i) EightZ2(i);...   % 7
    EightX1(i) EightY2(i) EightZ2(i)];     % 8

    EightVert = [EightVert; new];

    f = EightFaceBase+(8*(i-1));
    EightFace = [EightFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(5,2,8);
patch('Faces',EightFace,'Vertices',EightVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Eightboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.125]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Nine -----------------------
Nineboxestitle = 'k = 234';

boxesFileName = 'BivGaussian9.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Ninebase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

%Z1 = dataR(:,3)./dataR(:,2);  
NineZ1 = Ninebase;
NineZ2 = dataR(:,2)./dataR(:,1);% % height as count_in_box/volume
NineZ2 = NineZ2/CountsTotal; % height as relative_count_in_box/volume


NineX1 = dataR(:,3);
NineX2 = dataR(:,4);

NineY1 = dataR(:,5);
NineY2 = dataR(:,6);

Nineboxes = size(NineX1,1);

NineVert=[NineX1(1) NineY1(1) NineZ1(1);... % 1
    NineX2(1) NineY1(1) NineZ1(1);...   % 2
    NineX2(1) NineY2(1) NineZ1(1);...   % 3
    NineX1(1) NineY2(1) NineZ1(1);...   % 4
    NineX1(1) NineY1(1) NineZ2(1);...   % 5
    NineX2(1) NineY1(1) NineZ2(1);...   % 6
    NineX2(1) NineY2(1) NineZ2(1);...   % 7
    NineX1(1) NineY2(1) NineZ2(1)];     % 8

NineFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
NineFace = NineFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Nineboxes
    new = [NineX1(i) NineY1(i) NineZ1(i);... % 1
    NineX2(i) NineY1(i) NineZ1(i);...   % 2
    NineX2(i) NineY2(i) NineZ1(i);...   % 3
    NineX1(i) NineY2(i) NineZ1(i);...   % 4
    NineX1(i) NineY1(i) NineZ2(i);...   % 5
    NineX2(i) NineY1(i) NineZ2(i);...   % 6
    NineX2(i) NineY2(i) NineZ2(i);...   % 7
    NineX1(i) NineY2(i) NineZ2(i)];     % 8

    NineVert = [NineVert; new];

    f = NineFaceBase+(8*(i-1));
    NineFace = [NineFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(5,2,9);
patch('Faces',NineFace,'Vertices',NineVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Nineboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.125]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;

% --------------- histogram number Ten -----------------------
Tenboxestitle = 'k = 260';

boxesFileName = 'BivGaussian10.txt';

dataR = dlmread(boxesFileName, '\t', 0, 1); % miss out the labels

Tenbase = zeros(1,size(dataR,1));

Counts = dataR(:,2);
CountsTotal = sum(Counts);

%Z1 = dataR(:,3)./dataR(:,2);  
TenZ1 = Tenbase;
TenZ2 = dataR(:,2)./dataR(:,1); % height as count_in_box/volume
TenZ2 = TenZ2/CountsTotal; % height as relative_count_in_box/volume

TenX1 = dataR(:,3);
TenX2 = dataR(:,4);

TenY1 = dataR(:,5);
TenY2 = dataR(:,6);

Tenboxes = size(TenX1,1);

TenVert=[TenX1(1) TenY1(1) TenZ1(1);... % 1
    TenX2(1) TenY1(1) TenZ1(1);...   % 2
    TenX2(1) TenY2(1) TenZ1(1);...   % 3
    TenX1(1) TenY2(1) TenZ1(1);...   % 4
    TenX1(1) TenY1(1) TenZ2(1);...   % 5
    TenX2(1) TenY1(1) TenZ2(1);...   % 6
    TenX2(1) TenY2(1) TenZ2(1);...   % 7
    TenX1(1) TenY2(1) TenZ2(1)];     % 8

TenFaceBase = [1 2 3 4;...     % bottom
    5 6 7 8;...             % top
    1 2 6 5;...             % front
    3 4 8 7;...             % back
    2 3 7 6;...             % rhs    
    4 1 5 8];               % lhs
TenFace = TenFaceBase;

tcolbase = [1 0.6 0.7;...
    1 0.6 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7;...
    1 1 0.7];
tcol = tcolbase;

for i=2:Tenboxes
    new = [TenX1(i) TenY1(i) TenZ1(i);... % 1
    TenX2(i) TenY1(i) TenZ1(i);...   % 2
    TenX2(i) TenY2(i) TenZ1(i);...   % 3
    TenX1(i) TenY2(i) TenZ1(i);...   % 4
    TenX1(i) TenY1(i) TenZ2(i);...   % 5
    TenX2(i) TenY1(i) TenZ2(i);...   % 6
    TenX2(i) TenY2(i) TenZ2(i);...   % 7
    TenX1(i) TenY2(i) TenZ2(i)];     % 8

    TenVert = [TenVert; new];

    f = TenFaceBase+(8*(i-1));
    TenFace = [TenFace; f];

    t=tcolbase;
    tcol = [tcol; t];

end


subplot(5,2,10);
patch('Faces',TenFace,'Vertices',TenVert,'FaceVertexCData',tcol,...
      'FaceColor','flat','EdgeColor', color_edge)
alpha(0.1) % transparency
view(41.0, 34.0)
set(gca,'ZGrid','on')
title(Tenboxestitle, 'FontSize', 8, 'FontName', 'Ariel', 'FontWeight', 'Bold')   
zlim([0 0.125]); 
camlight right; lighting phong;
% PDF
hold on;
surface(SGx,SGy,SGz,'FaceColor','blue','EdgeColor','none');
alpha(.5) % transparency
camlight right; lighting phong;
