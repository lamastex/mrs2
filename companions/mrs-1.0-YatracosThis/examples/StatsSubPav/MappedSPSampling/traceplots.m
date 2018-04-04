close all
clear all
clc

tol=0.1;

figure
subplot(2,1,1)
GRDiag = importdata(['GelmanRubinLeavesScalar.txt']);
GRDiag = GRDiag.data;

Converged=55510;
fs=17;

%% Convergence at ...
rhatFlag = GRDiag(:,5);
%rhatSampled = GRDiag(:,7);


loglog([Converged Converged], [min(GRDiag(:,5)+0.001) max(GRDiag(:,5))], 'k--')
hold on
loglog([1 length(GRDiag)], [1-tol 1-tol], 'k-')
loglog([1 length(GRDiag)], [1+tol 1+tol], 'k-')
loglog(GRDiag(:,1),GRDiag(:,5), 'k-.', 'MarkerFaceColor', 'k', 'LineWidth', 2)
%axis tight
xlabel('Time step', 'FontSize', fs);
ylabel('Rhat', 'FontSize', fs)
%title('Gelman-Rubin Convergence Diagnostics for #Leaves')
get(0, 'CurrentFigure')
set(gca, 'FontSize', fs)



%% Trace loglogs of scalar summary
subplot(2,1,2)
fileName = ['LeavesScalar.txt'];
leafSummary = importdata(fileName);
leafSummary = leafSummary.data;
semilogx(1:length(leafSummary), leafSummary(:,2), 'k--', 'LineWidth', 2)
hold on
semilogx(1:length(leafSummary), leafSummary(:,3), 'k') %, 'LineWidth', 2)
semilogx([Converged Converged], [-200+min(min(leafSummary(:,2), leafSummary(:,3))) max(max(leafSummary(:,2)), max(leafSummary(:,3)))], 'k--')
axis([1 length(leafSummary) -200 200+(max(leafSummary(:,3)))])


xlabel('Time step', 'FontSize', fs)
ylabel('#leaves', 'FontSize', fs)
%legend('Chain 1', 'Chain 2', 'Location', 'SouthEast')
get(0, 'CurrentFigure')
set(gca, 'FontSize', fs)
print(gcf, '-depsc', ...
['~/svn_mrs/mrs/branches/gat41/manuscript/Figures/Gaussian2DGRTraceplots.eps']);


M = max(leafSummary)

find(leafSummary(:,3)==M(3))

figure
myhist = importdata('EstFunctionAfterNormalized.txt');
myhist = myhist.data;
plotEstFunction2D(myhist, 0, 1)
get(0, 'CurrentFigure')
set(gca, 'FontSize', fs)
axis tight
print(gcf, '-depsc', ...
['~/svn_mrs/mrs/branches/gat41/manuscript/Figures/Gaussian2DEstHist.eps']);



figure
myhist = importdata('FinalHist0.txt');
myhist = myhist.data;
plotEstFunction2D(myhist, 0, 1)
get(0, 'CurrentFigure')
set(gca, 'FontSize', fs)
axis tight
print(gcf, '-depsc', ...
['~/svn_mrs/mrs/branches/gat41/manuscript/Figures/Gaussian2DGRMCMCHist.eps']);
