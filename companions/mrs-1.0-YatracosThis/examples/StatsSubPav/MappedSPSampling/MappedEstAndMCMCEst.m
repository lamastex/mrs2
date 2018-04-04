%% M-file to plot the traceplots, and Gelman-Rubin convergence diagnostics
%% for Gaussian.
clc
clear all
close all
 
%density1 = 'Rosen2D';
%density2 = 'Rosen2Da';

density1 = 'Gaussian2D';
density2 = 'Gaussian2Da';

% density1 = 'Gaussian1D';
% density2 = 'Gaussian1Da';

%D=1;
D=2;

fontsize=17.5;
fontsize2 = 17.5;

% user inputs
tol = 0.1;

%% plot histograms
figure('OuterPosition', [900 900 900 900])

%estimated function
%myHist = importdata('EstFunction.txt');
myHist = importdata(['NewDataForPaper/EstFunction', density1, '.txt']);
myHist = myHist.data;
totalArea = 1;

if ( D == 1)
    plotEstFunction1D(myHist, 2, 'k', 1);
    xlabel('x', 'FontSize', fontsize)
    ylabel('Density', 'FontSize', fontsize)
    
else
    plotEstFunction2D(myHist, 0, totalArea);
    xlabel('x_1', 'FontSize', fontsize)
    ylabel('x_2', 'FontSize', fontsize)
    zlabel('Density', 'FontSize', fontsize)

    get(0, 'CurrentFigure')
    set(gca, 'FontSize', fontsize)

    %gat41@ubuntu
    saveas(gcf, ...
    ['~/svn_mrs/mrs/branches/gat41/manuscript/Figures/newMCMCHistogramActual', ...
    density1, '.eps']);
end

avgHist = importdata(['NewDataForPaper/FinalHist', density1, '.txt']);
%avgHist = importdata('FinalHist.txt');
avgHist = avgHist.data;
totalArea = 1;

if (D==1)
    hold on 
    plotEstFunction1D(avgHist, 2, 'k', 2);
    xlabel('x', 'FontSize', fontsize)
    ylabel('Density estimate', 'FontSize', fontsize)
    
else
    figure('OuterPosition', [900 900 900 900])
    plotEstFunction2D(avgHist, 0, 1);
    xlabel('x_1', 'FontSize', fontsize)
    ylabel('x_2', 'FontSize', fontsize)
    zlabel('Density estimate', 'FontSize', fontsize)
    
    get(0, 'CurrentFigure')
    set(gca, 'FontSize', fontsize)

    %gat41@ubuntu
    saveas(gcf, ...
    ['~/svn_mrs/mrs/branches/gat41/manuscript/Figures/newMCMCHistogram1', ...
    density1, '.eps']);    
end


avgHist = importdata(['NewDataForPaper/FinalHist', density2, '.txt']);
%avgHist = importdata('FinalHist.txt');
avgHist = avgHist.data;
totalArea = 1;

if (D==1)
    hold on
    plotEstFunction1D(avgHist, 2, [0.5 0.5 0.5], 2);
    xlabel('x', 'FontSize', fontsize)
    ylabel('Density estimate', 'FontSize', fontsize)
    get(0, 'CurrentFigure')
    set(gca, 'FontSize', fontsize)

%gat41@ubuntu
saveas(gcf, ...
    ['~/svn_mrs/mrs/branches/gat41/manuscript/Figures/newMCMCHistogram', ...
    density1, '.eps']);    

else
    figure('OuterPosition', [900 900 900 900])
    plotEstFunction2D(avgHist, 0, 1);
    xlabel('x_1', 'FontSize', fontsize)
    ylabel('x_2', 'FontSize', fontsize)
    zlabel('Density estimate', 'FontSize', fontsize)
    
    get(0, 'CurrentFigure')
    set(gca, 'FontSize', fontsize)

    %gat41@ubuntu
    saveas(gcf, ...
    ['~/svn_mrs/mrs/branches/gat41/manuscript/Figures/newMCMCHistogram2', ...
    density1, '.eps']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trace plots
%density 1 
%figure('OuterPosition', [1600 2100 1600 1100])
figure 
subplot(2,1,2)
GRDiag = importdata(['NewDataForPaper/Gelman', density1, '.txt']);
%GRDiag = importdata('GelmanRubinLeavesScalar.txt');
GRDiag = GRDiag.data;
loglog(GRDiag(:,1),GRDiag(:,5), 'k-.', 'MarkerFaceColor', 'k', 'LineWidth', 2)

% Convergence at ...
rhatFlag = GRDiag(:,6);
rhatSampled = GRDiag(:,7);
Converged1 = min(find (rhatFlag==1 & rhatSampled==1));
R1 = GRDiag(Converged1, 5); 


L1 = length(GRDiag);
M1 = max(GRDiag(:,5));

% %density 2
GRDiag = importdata(['NewDataForPaper/Gelman', density2, '.txt']);
%GRDiag = importdata('GelmanRubinLeavesScalar.txt');
GRDiag = GRDiag.data;
hold on
loglog(GRDiag(:,1),GRDiag(:,5), '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 2)

rhatFlag = GRDiag(:,6);
rhatSampled = GRDiag(:,7);
Converged2 = min(find (rhatFlag==1 & rhatSampled==1));
R2 = GRDiag(Converged2, 5);


L2 = length(GRDiag);
M2 = max(GRDiag(:,5));

legend('H1', 'H2')

L = max(L1, L2);
M = max(M1, M2);

loglog([1 L], [1-tol 1-tol], 'k-')
hold on
loglog([1 L], [1+tol 1+tol], 'k-')
loglog(Converged1, R1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize',7)
hold on
loglog(Converged2, R2, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 7)

axis([1 L+10000, 1-tol-0.5 M+1000])
set(gca, 'YTick', logspace(1,7,7))

xlabel('Time step', 'FontSize', fontsize2);
ylabel('Rhat', 'FontSize', fontsize2)
get(0, 'CurrentFigure')
set(gca, 'FontSize', fontsize2)

%% Trace plots of scalar summary
subplot(2,1,1)
%density 1
fileName = ['NewDataForPaper/LeavesScalar', density1, '.txt'];
%fileName = 'LeavesScalar.txt';
leafSummary = importdata(fileName);
leafSummary = leafSummary.data;
loglog(1:length(leafSummary), leafSummary(:,2), 'k-x', 'LineWidth', 2)
hold on
loglog(1:length(leafSummary), leafSummary(:,3), 'k-', 'LineWidth', 2)

M1 = max(leafSummary(:,3));


%density 2
fileName = ['NewDataForPaper/LeavesScalar', density2, '.txt'];
%fileName = 'LeavesScalar.txt';
leafSummary = importdata(fileName);
leafSummary = leafSummary.data;
loglog(1:length(leafSummary), leafSummary(:,2), '--', 'LineWidth', 2, ...
    'Color', [0.5 0.5 0.5])
hold on
loglog(1:length(leafSummary), leafSummary(:,3), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2)

M2 = max(leafSummary(:,3));

M = max(M1, M2);

%loglog(Converged1, R1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3)
%hold on
%loglog(Converged2, R2, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], ...
%    'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 3)

%axis auto
axis([1 L+10000, 1 M+10000])

xlabel('Time step', 'FontSize', fontsize2)
ylabel('Number of leaves', 'FontSize', fontsize2)
%legend('C1 (H1)', 'C2 (H1)', 'C1 (H2)', 'C2 (H2)', 'Location', 'NorthEast')

get(0, 'CurrentFigure')

set(gca, 'YTick', logspace(1,10,10))

set(gca, 'FontSize', fontsize2)
print(gcf, '-depsc', ...
['~/svn_mrs/mrs/branches/gat41/manuscript/Figures/newGRPlots', density1, '.eps']);
