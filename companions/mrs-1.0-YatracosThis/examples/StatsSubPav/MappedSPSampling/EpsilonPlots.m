clear all
close all
clc

figure
eps = importdata('Eps.txt');
plot(1:1:length(eps), eps, 'LineWidth', 2)
hold on
% loglog(10, eps(10), 'bo')
% loglog(100, eps(100), 'bo')
% loglog(100, eps(1000), 'bo')
% loglog(100, eps(1000), 'bo')

% eps = importdata('Epsilon/EpsGaussian2D.txt');
% loglog(1:1:length(eps), eps, 'k', 'LineWidth', 2)
% % 
% 
% 
% %loglog([1 length(eps)], [0.5 0.5])
% 
% % 
% eps = importdata('Epsilon/EpsGaussian10D.txt');
% loglog(1:1:length(eps), eps, 'Color', [0.5 0.5 0.5], 'LineWidth', 2)
% % 
% 
% % 
% %%
% 
% eps = importdata('Epsilon/EpsRosen2D.txt');
% loglog(1:1:length(eps), eps, 'g', 'LineWidth', 2)
% 
% hold on
% 
% eps = importdata('Epsilon/EpsRosen10D.txt');
% loglog(1:1:length(eps), eps, 'r', 'LineWidth', 2)
% 
% 
% legend('Location', 'SouthWest', 'G1D', 'G2D', 'G10D', 'R2D', 'R10D')
% 
% xlabel('\Lambda', 'FontSize', 16)
% ylabel('\epsilon', 'FontSize', 16)
% % 
% get(0, 'CurrentFigure')
% set(gca, 'FontSize', 16)
% 
% %get(0, 'CurrentFigure')
% %print(gcf, '-dpng', ...
% %['~/svn_mrs/mrs/branches/gat41/manuscript/', density, '.png']);
% 
