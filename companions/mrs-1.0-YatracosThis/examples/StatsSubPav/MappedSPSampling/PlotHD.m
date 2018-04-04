clc 
clear all
close all

%Plot Hellinger distances (avg and standard deviation)

L = [10^2, 10^3, 10^4, 10^5, 10^6];
mu = zeros(5,1);

G1cilow = zeros(5,1);
G1ciupp = zeros(5,1);

G2cilow = zeros(5,1);
G2ciupp = zeros(5,1);

G5cilow = zeros(5,1);
G5ciupp = zeros(5,1);

R1cilow = zeros(5,1);
R1ciupp = zeros(5,1);

R2cilow = zeros(5,1);
R2ciupp = zeros(5,1);

R5cilow = zeros(5,1);
R5ciupp = zeros(5,1);


%% Gaussian1D
%load Gaussian2DHD/Gaussian2DL100Actual
load Gaussian1DL100Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(1) = muh;
G1cilow(1) = muci(1);
G1ciupp(1) = muci(2);

load Gaussian1DL1000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(2) = muh;
G1cilow(2) = muci(1);
G1ciupp(2) = muci(2);

load Gaussian1DL10000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(3) = muh;
G1cilow(3) = muci(1);
G1ciupp(3) = muci(2);

load Gaussian1DL100000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(4) = muh;
G1cilow(4) = muci(1);
G1ciupp(4) = muci(2);

load Gaussian1DL1000000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(5) = muh;
G1cilow(5) = muci(1);
G1ciupp(5) = muci(2);

loglog(L, mu, '-ko', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 2);
hold on


%% Gaussian2D
load Gaussian2DL100Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(1) = muh;
G2cilow(1) = muci(1);
G2ciupp(1) = muci(2);

load Gaussian2DL1000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(2) = muh;
G2cilow(2) = muci(1);
G2ciupp(2) = muci(2);

load Gaussian2DL10000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(3) = muh;
G2cilow(3) = muci(1);
G2ciupp(3) = muci(2);

load Gaussian2DL100000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(4) = muh;
G2cilow(4) = muci(1);
G2ciupp(4) = muci(2);

load Gaussian2DL1000000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(5) = muh;
G2cilow(5) = muci(1);
G2ciupp(5) = muci(2);

%% Gaussian5D
load Gaussian5DL100Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(1) = muh;
G5cilow(1) = muci(1);
G5ciupp(1) = muci(2);

load Gaussian5DL1000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(2) = muh;
G5cilow(2) = muci(1);
G5ciupp(2) = muci(2);

load Gaussian5DL10000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(3) = muh;
G5cilow(3) = muci(1);
G5ciupp(3) = muci(2);

load Gaussian5DL100000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(4) = muh;
G5cilow(4) = muci(1);
G5ciupp(4) = muci(2);

load Gaussian5DL1000000Actual
[muh sigma muci sigmaci] = normfit(HDActual);
mu(5) = muh;
G5cilow(5) = muci(1);
G5ciupp(5) = muci(2);

loglog(L, mu, '-ko', 'LineWidth', 3, 'MarkerFaceColor', 'k', 'MarkerSize', 2);


%% Rosen1D
load Rosen1DL100TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(1) = muh;
R1cilow(1) = muci(1);
R1ciupp(1) = muci(2);

load Rosen1DL1000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(2) = muh;
R1cilow(2) = muci(1);
R1ciupp(2) = muci(2);

load Rosen1DL10000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(3) = muh;
R1cilow(3) = muci(1);
R1ciupp(3) = muci(2);

load Rosen1DL100000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(4) = muh;
R1cilow(4) = muci(1);
R1ciupp(4) = muci(2);

load Rosen1DL1000000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(5) = muh;
R1cilow(5) = muci(1);
R1ciupp(5) = muci(2);

loglog(L, mu, '--ko', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 2);

%% Rosen2D
load Rosen2DL100TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(1) = muh;
R2cilow(1) = muci(1);
R2ciupp(1) = muci(2);

load Rosen2DL1000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(2) = muh;
R2cilow(2) = muci(1);
R2ciupp(2) = muci(2);

load Rosen2DL10000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(3) = muh;
R2cilow(3) = muci(1);
R2ciupp(3) = muci(2);

load Rosen2DL100000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(4) = muh;
R2cilow(4) = muci(1);
R2ciupp(4) = muci(2);

load Rosen2DL1000000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(5) = muh;
R2cilow(5) = muci(1);
R2ciupp(5) = muci(2);

loglog(L, mu, '--ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 2);


%% Rosen5D
load Rosen5DL100TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(1) = muh;
R5cilow(1) = muci(1);
R5ciupp(1) = muci(2);

load Rosen5DL1000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(2) = muh;
R5cilow(2) = muci(1);
R5ciupp(2) = muci(2);

load Rosen5DL10000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(3) = muh;
R5cilow(3) = muci(1);
R5ciupp(3) = muci(2);

load Rosen5DL100000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(4) = muh;
R5cilow(4) = muci(1);
R5ciupp(4) = muci(2);

load Rosen5DL1000000TrueActual
[muh sigma muci sigmaci] = normfit(HDTrueActual);
mu(5) = muh;
R5cilow(5) = muci(1);
R5ciupp(5) = muci(2);

legend('1D-Gauss', '2D-Gauss', '5D-Gauss', '1D-Rosen', '2D-Rosen', '5D-Rosen')


%% plot the confidence intervals
loglog(L, G1cilow, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
loglog(L, G1ciupp, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5],'MarkerEdgeColor', [0.5, 0.5, 0.5]);

for i = 1:length(L)
    loglog([L(i) L(i)], [G1cilow(i) G1ciupp(i)], 'k')
end

loglog(L, mu, '-ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 2);

loglog(L, G2cilow, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
loglog(L, G2ciupp, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5],'MarkerEdgeColor', [0.5, 0.5, 0.5]);

for i = 1:length(L)
    loglog([L(i) L(i)], [G2cilow(i) G2ciupp(i)], 'k')
end

loglog(L, G5cilow, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
loglog(L, G5ciupp, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5],'MarkerEdgeColor', [0.5, 0.5, 0.5]);

for i = 1:length(L)
    loglog([L(i) L(i)], [G5cilow(i) G5ciupp(i)], 'k')
end

loglog(L, R1cilow, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
loglog(L, R1ciupp, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5],'MarkerEdgeColor', [0.5, 0.5, 0.5]);

for i = 1:length(L)
    loglog([L(i) L(i)], [R1cilow(i) R1ciupp(i)], 'k')
end

loglog(L, R2cilow, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
loglog(L, R2ciupp, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5],'MarkerEdgeColor', [0.5, 0.5, 0.5]);

for i = 1:length(L)
    loglog([L(i) L(i)], [R2cilow(i) R2ciupp(i)], 'k')
end

loglog(L, mu, '--ko', 'LineWidth', 3, 'MarkerFaceColor', 'k', 'MarkerSize', 2);

loglog(L, R5cilow, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
loglog(L, R5ciupp, 'o', 'MarkerSize', 2', 'MarkerFaceColor', [0.5, 0.5, 0.5],'MarkerEdgeColor', [0.5, 0.5, 0.5]);

for i = 1:length(L)
    loglog([L(i) L(i)], [R5cilow(i) R5ciupp(i)], 'k')
end

get(0, 'CurrentFigure')
set(gca, 'FontSize', 16)
xlabel('\Lambda', 'FontSize', 16)
ylabel('Hellinger distance', 'FontSize', 16)

% print(gcf, '-depsc', ...
%  ['mrs/branches/gat41/manuscript/Figures/HDCIPlots.eps']);
