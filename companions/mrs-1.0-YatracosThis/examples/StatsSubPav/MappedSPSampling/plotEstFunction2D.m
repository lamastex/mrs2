function [fhat xlow xupp ylow yupp]=plotEstFunction2D(data, hz,totalArea)
%           
%Purpose: To plot the histogram using subpaving boxes.
%
%Input argument list:
%       data:   data set
%       N:      total number of points
%       hz  :   position of subpaving plot
%
%Output argument list:
%       fhat: density estimate
%       Lb  : lower bound
%       Ub  : upper bound
%       vol  : the Lebesgue measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volBin=data(:,1);
xlow=data(:,3)';
xupp=data(:,4)';
ylow=data(:,5)';
yupp=data(:,6)';
fhat=data(:,2)/totalArea;
fhat=fhat';


%base and fhat (top and bottom)
XX=[1 0 0 1]'*xlow + [0 1 1 0]'*xupp; 
YY=[1 1 0 0 ]'*ylow + [0 0 1 1]'*yupp;
ZZ=[zeros(4,1)]*fhat;
fill3(XX, YY, ZZ, 'w','FaceAlpha', .25)
%fill3(XX, YY, ZZ, 'k','FaceAlpha', .25)

grid on
%axis([min(xlow) max(xupp) min(ylow) max(yupp) hz max(fhat)])
hold on
XX=[1 0 0 1]'*xlow + [0 1 1 0]'*xupp; 
YY=[1 1 0 0 ]'*ylow + [0 0 1 1]'*yupp;
ZZ=[ones(4,1)]*fhat;
fill3(XX, YY, ZZ, 'w','FaceAlpha', .25)
%fill3(XX, YY, ZZ, 'k','FaceAlpha', .25)

%front and back face
XX=[ones(4,1)]*xlow + [zeros(4,1)]*xupp;
YY=[0 1 1 0]'*ylow + [1 0 0 1]'*yupp;
ZZ=[0 0 1 1]'*fhat;
fill3(XX, YY, ZZ, 'w','FaceAlpha', .25)
%fill3(XX, YY, ZZ, 'k','FaceAlpha', .25)

XX=[zeros(4,1)]*xlow + [ones(4,1)]*xupp;
YY=[0 1 1 0]'*ylow + [1 0 0 1]'*yupp;
ZZ=[0 0 1 1]'*fhat;
fill3(XX, YY, ZZ, 'w','FaceAlpha', .25)
%fill3(XX, YY, ZZ, 'k','FaceAlpha', .25)

%left and right face
YY=[ones(4,1)]*ylow + [zeros(4,1)]*yupp;
XX=[0 1 1 0]'*xlow + [1 0 0 1]'*xupp;
ZZ=[0 0 1 1]'*fhat;
fill3(XX, YY, ZZ, 'w','FaceAlpha', .25)
%fill3(XX, YY, ZZ, 'k','FaceAlpha', .25)

YY=[zeros(4,1)]*ylow + [ones(4,1)]*yupp;
XX=[0 1 1 0]'*xlow + [1 0 0 1]'*xupp;
ZZ=[0 0 1 1]'*fhat;
fill3(XX, YY, ZZ, 'w','FaceAlpha', .25)
%fill3(XX, YY, ZZ, 'k','FaceAlpha', .25)

xlabel('x', 'FontSize', 18)
ylabel('y', 'FontSize', 18)
zlabel('f_{n, s}', 'FontSize', 18)
get(0, 'CurrentFigure')
set(gca, 'FontSize', 18)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the default 3-D perspective view
%view(3)
view(-36.0, 25.0)
%End of function: subpavhist2D.m