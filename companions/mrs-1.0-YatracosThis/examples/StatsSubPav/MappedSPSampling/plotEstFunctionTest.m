clc
clear all
close all

%1D
subplot(1,2,1)
myHist = importdata('EstFunction.txt');
%myHist = importdata('EstFunctionAfterNormalized.txt');
myHist = myHist.data;
totalArea = 1;
%plotEstFunction1D(myHist, totalArea);
plotEstFunction2D(myHist, 0, totalArea);
title('MappedSP Function')

subplot(1,2,2)
%myHist = importdata('EstFunction.txt');
myHist = importdata('EstFunctionAfterNormalized.txt');
myHist = myHist.data;
totalArea = 1;
%plotEstFunction1D(myHist, totalArea);
plotEstFunction2D(myHist, 0, totalArea);
title('MappedSP Function')


% subplot(1,2,2)
% avgHist = importdata('FinalHist.txt');
% avgHist = avgHist.data;
% totalArea = 1;
% %plotEstFunction1D(avgHist, 1);
% plotEstFunction2D(avgHist, 0, 1);
% title('MCMC Averaged Hist')

% testHist = importdata('HistTest.txt');
% testHist = testHist.data;
% plotEstFunction1D(testHist);
