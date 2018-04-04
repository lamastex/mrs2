%% plots for PriorityQueue IAE (FinMix)
clear all
close all
clc
addpath ../


startN = 1;
endN = 1;

a=1; %n=10^4
b=0; %n=10^5
c=0; %n=10^6

plotting=1;

numhist=1;

mix = 1;


if (a==1) count1 = 1:1:100; end
if (b==1) count2 = 1:1:100; end
if (c==1) count3 = 1:1:100000; end


%% container to keep optimal bin size
minMIAEPQ=[];
minSDPQ=[];
minBinPQ=[];

minMIAEMDE=[];
minSDMDE=[];
minBinMDE=[];

minMIAEReg=[];
minSDReg=[];
minBinReg=[];


if (a==1)
    %============================n=10^4=====================================%
    load BinsPQ10000
    load BinsReg10000
    load IAEPQ10000
    load IAEReg10000
    
    %optimal bin size for PQ
    MIAEPQ = mean(IAEPQ);
    minMIAEPQ(1) = min(MIAEPQ);
    binIndex= find(MIAEPQ==minMIAEPQ(1));
    minBinPQ = binsPQ(:,min(binIndex)); %would be a range  
    minSDPQ(1) = std(IAEPQ(:,min(binIndex)));
    
    %is this a fair comparison?
    %optimal bin size for reg
    minMIAEReg(1) = mean(IAEReg);
    minSDReg(1) = std(IAEReg);
    minBinReg(1) = min(binsReg);
    
    %% plots for PQ
    figure
    subplot(3,1,1)
    semilogy(binsPQ, IAEPQ, 'b.')
    title('IAE vs number of bins')
    xlabel('bins')
    ylabel('IAE')
    
    subplot(3,1,2)
    plot(count1, IAEPQ, 'b.')
    hold on
    plot(count1, MIAEPQ, 'r.')
    title('IAE vs counts')
    xlabel('counts')
    ylabel('IAE')
    
    subplot(3,1,3)
    plot(binsReg, IAEReg, 'b.')
    title('RegHist: IAE vs num of bins')
    xlabel('number of bins')
    ylabel('IAE')    
    
    %% MDE All
    minTheta = [];
    minIAEBin = [];
    thetaIndexCount = 1;
    for h = 1:numhist
        %read the data
        myFile=['DeltaMax', int2str(h), '.txt'];
        delta = dlmread(myFile);
        [nrow ncol] = size(delta);
                
        %plot delta for hist h
        if (plotting == 1)
            figure
        end
        
        %plot IAE
        IAEMDE = importdata(['IAEHistMDE', int2str(h), '.txt']);
        minIAEBin(h) = min(find(IAEMDE==min(IAEMDE)));
        if (plotting == 1)
            title(['\Delta_\theta vs theta for ', myFile])
            xlabel('theta')
            ylabel('\Delta_\theta (linear scale) (blue) /IAE (green)')
            hold on
            plot(0:1:(length(IAEMDE)-1), IAEMDE, '-gd', 'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'g')
        end
        
        % plot delta max
        seeInf=[];
        for i=2:nrow
            theta = 0:1:(i-1);
            dd = delta(i, 1:i);
            infDD=min(dd);
            infTheta=theta(find(dd==infDD));
            seeInf(i-1) = min(infTheta);
            
            if (plotting==1)
                %plot delta
                plot(theta, dd, '-bx', 'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 10)
                hold on
                %plot the infimum
                plot(infTheta, infDD, 'o', 'MarkerFaceColor', 'r', ...
                    'MarkerEdgeColor', 'k')
             end
        end
          
        %get the thetas that gives the minimum delta at the 'final' stage
        thetaStar(h) = min(infTheta);
%         figure
%         plot(1:(nrow-1), seeInf, 'bo', 'MarkerFaceColor', 'b')
        [u trash loc] = unique(seeInf);
        freq = histc(loc, 1:length(u));
        maxFreqIndex = find(freq==max(freq));
        allInf(h) = min(u(maxFreqIndex));
%         hold on
%         plot(1, u(maxFreqIndex), 'ro', 'MarkerFaceColor', 'r')
%         
    end
    
    %comparison of optimal bin
    figure
    plot([0 100], [0 100])
    hold on
    plot(thetaStar, minIAEBin, 'bo', 'MarkerFaceColor', 'b');
    title('The Optimal Bin')
    xlabel('Bin from MDE')
    ylabel('Bin from IAE')
    
    %frequency count for thetaStar
    uniqMinTheta = unique(thetaStar);
    uniqIAE = unique(minIAEBin);
    freqMDE=[];
    freqIAE=[];
    for u = 1:length(uniqMinTheta)
        freqMDE(u) = sum((thetaStar == uniqMinTheta(u)));
    end
    for u = 1:length(uniqIAE)
        freqIAE(u) = sum(minIAEBin==uniqIAE(u));
    end

    minBinMDE(1) = uniqMinTheta(min(find(freqMDE == max(freqMDE))))+1;
    hold on
    plot(minBinMDE(1), uniqIAE(min(find(freqIAE == max(freqIAE)))), 'ro', 'MarkerFaceColor', 'r');
    
   [u trash loc] = unique(allInf);
   freq = histc(loc, 1:length(u));
     
    
    
    [uniqMinTheta+1; freqMDE]
    [u+1; freq] 
    [uniqIAE; freqIAE]

    %MIAE and SD for this bin
    %ind = uniqMinTheta(min(find(freq == max(freq))))+1;
    
    %minMIAEMDE(1)=mean(IAEminTheta);
    %minSDMDE (1)=std(IAEminTheta);
    %============================n=10^4=====================================%
end

if (b==1)
    %============================n=10^5=====================================%
    load BinsPQ100000
    load BinsReg100000
    load IAEPQ100000
    load IAEReg100000
    
    
    %optimal bin size for PQ
    MIAEPQ = mean(IAEPQ);
    minMIAEPQ(2) = min(MIAEPQ);
    binIndex= find(MIAEPQ==minMIAEPQ(2));
    minBinPQ = binsPQ(:,min(binIndex)); %would be a range  
    minSDPQ(2) = std(IAEPQ(:,min(binIndex)));
    
    %is this a fair comparison?
    %optimal bin size for reg
    minMIAEReg(2) = mean(IAEReg);
    minSDReg(2) = std(IAEReg);
    minBinReg(2) = min(binsReg);
    
    %% plots for PQ
    figure
    subplot(3,1,1)
    semilogy(binsPQ, IAEPQ, 'b.')
    title('IAE vs number of bins')
    xlabel('bins')
    ylabel('IAE')
    
    subplot(3,1,2)
    plot(count2, IAEPQ, 'b.')
    hold on
    plot(count2, MIAEPQ, 'r.')
    title('IAE vs counts')
    xlabel('counts')
    ylabel('IAE')
    
    subplot(3,1,3)
    plot(binsReg, IAEReg, 'b.')
    title('RegHist: IAE vs num of bins')
    xlabel('number of bins')
    ylabel('IAE')    
    
    %% MDE All
    minTheta = [];
    minIAEBin = [];
    thetaIndexCount = 1;
    for h = 1:numhist
        %read the data
        myFile=['n5/DeltaMax', int2str(h), '.txt'];
        delta = dlmread(myFile);
        [nrow ncol] = size(delta);
                
        %plot delta for hist h
        if (plotting == 1)
            figure
        end
        
        %plot IAE
        IAEMDE = importdata(['n5/IAEHistMDE', int2str(h), '.txt']);
        minIAEBin(h) = min(find(IAEMDE==min(IAEMDE)));
        if (plotting == 1)
            title(['\Delta_\theta vs theta for ', myFile])
            xlabel('theta')
            ylabel('\Delta_\theta (linear scale) (blue) /IAE (green)')
            hold on
            plot(0:1:(length(IAEMDE)-1), IAEMDE, '-gd', 'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'g')
        end
        
        % plot delta max
        seeInf=[];
        for i=2:nrow
            theta = 0:1:(i-1);
            dd = delta(i, 1:i);
            infDD=min(dd);
            infTheta=theta(find(dd==infDD));
            seeInf(i-1) = min(infTheta);
            
            if (plotting==1)
                %plot delta
                plot(theta, dd, '-bx', 'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 10)
                hold on
                %plot the infimum
                plot(infTheta, infDD, 'o', 'MarkerFaceColor', 'r', ...
                    'MarkerEdgeColor', 'k')
             end
        end
          
        %get the thetas that gives the minimum delta at the 'final' stage
        thetaStar(h) = min(infTheta);
%         figure
%         plot(1:(nrow-1), seeInf, 'bo', 'MarkerFaceColor', 'b')
        [u trash loc] = unique(seeInf);
        freq = histc(loc, 1:length(u));
        maxFreqIndex = find(freq==max(freq));
        allInf(h) = min(u(maxFreqIndex));
%         hold on
%         plot(1, u(maxFreqIndex), 'ro', 'MarkerFaceColor', 'r')
%         
    end
    
    %comparison of optimal bin
    figure
    plot([0 100], [0 100])
    hold on
    plot(thetaStar, minIAEBin, 'bo', 'MarkerFaceColor', 'b');
    title('The Optimal Bin')
    xlabel('Bin from MDE')
    ylabel('Bin from IAE')
    
    %frequency count for thetaStar
    uniqMinTheta = unique(thetaStar);
    uniqIAE = unique(minIAEBin);
    freqMDE=[];
    freqIAE=[];
    for u = 1:length(uniqMinTheta)
        freqMDE(u) = sum((thetaStar == uniqMinTheta(u)));
    end
    for u = 1:length(uniqIAE)
        freqIAE(u) = sum(minIAEBin==uniqIAE(u));
    end

    minBinMDE(1) = uniqMinTheta(min(find(freqMDE == max(freqMDE))))+1;
    hold on
    plot(minBinMDE(1), uniqIAE(min(find(freqIAE == max(freqIAE)))), 'ro', 'MarkerFaceColor', 'r');
    
   [u trash loc] = unique(allInf);
   freq = histc(loc, 1:length(u));
     
    
    
    [uniqMinTheta+1; freqMDE]
    [u+1; freq] 
    [uniqIAE; freqIAE]

    %MIAE and SD for this bin
    %ind = uniqMinTheta(min(find(freq == max(freq))))+1;
    
    %minMIAEMDE(1)=mean(IAEminTheta);
    %minSDMDE (1)=std(IAEminTheta);
    %============================n=10^4=====================================%
end

if (c==1)
    %============================n=10^6=====================================%
    load BinsPQ1000000
    load BinsReg1000000
    load IAEPQ1000000
    load IAEReg1000000
    
    
    %optimal bin size for PQ
    MIAEPQ = mean(IAEPQ);
    minMIAEPQ(3) = min(mean(IAEPQ));
    binIndex= find(mean(IAEPQ)==minMIAEPQ(3));
    minBin = binsPQ(:,min(binIndex));
    minBinPQ(3) = floor(mean(minBin)); %a bit vary about this
    %minBinPQ(1) = min(minBin); %a bit vary about this
    minSDPQ(3) = std(IAEPQ(:,min(binIndex)));
    
    %is this a fair comparison?
    %optimal bin size for reg
    minMIAEReg(3) = mean(IAEReg);
    minSDReg(3) = std(IAEReg);
    minBinReg(3) = min(binsReg);
    
    %plots
    figure
    subplot(3,1,1)
    semilogy(binsPQ, IAEPQ, 'b.')
    title('IAE vs number of bins')
    xlabel('bins')
    ylabel('IAE')
    
    subplot(3,1,2)
    plot(count3, IAEPQ, 'b.')
    hold on
    plot(count3, MIAEPQ, 'r.')
    title('IAE vs counts')
    xlabel('counts')
    ylabel('IAE')
    
    subplot(3,1,3)
    plot(binsReg, IAEReg, 'b.')
    title('RegHist: IAE vs num of bins')
    xlabel('number of bins')
    ylabel('IAE')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MDE All
    minTheta = [];
    for h = 1:numhist
        myFile=['FM', int2str(mix), '/n6/DeltaMax', int2str(h), '.txt'];
        IAEMDE = importdata(['FM', int2str(mix), '/n6/IAEHistMDE', int2str(h), '.txt']);
        delta = dlmread(myFile);
        [nrow ncol] = size(delta);
        
        
        if (plotting == 1)
            figure
        end
        
        %plot delta
        for i=2:nrow
            theta = 0:1:(i-1);
            dd = delta(i, 1:i);
            infDD=min(dd);
            infTheta=theta(find(dd==infDD));
            
            if (plotting == 1)
                plot(theta, dd, '-bx', 'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 10)
                hold on
                %plot the infimum
                plot(infTheta, infDD, 'o', 'MarkerFaceColor', 'r', ...
                    'MarkerEdgeColor', 'k')
            end
        end
        
        if (plotting == 1)
            title(['\Delta_\theta vs theta for ', myFile])
            xlabel('theta')
            ylabel('\Delta_\theta (linear scale) (blue) /IAE (green)')
            
            %plot IAE
            hold on
            plot(0:1:(length(IAEMDE)-1), IAEMDE, '-gd', 'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'g')
        end
        %thetaStar(h)=min(theta(find(dd < (infDD + 1/10^4))));
        thetaStar(h)=min(infTheta);
        minIAEMDE(h)=min(find(IAEMDE==min(IAEMDE)));
    end
    
    uniqMinTheta = unique(thetaStar);
    freq=[];
    for u = 1:length(uniqMinTheta)
        freq(u) = sum((thetaStar == uniqMinTheta(u)));
    end
    minBinMDE(3) = uniqMinTheta(min(find(freq == max(freq))))+1;
    
    %MIAE and SD for this bin
    ind = uniqMinTheta(min(find(freq == max(freq))))+1;
    IAEminTheta = [];
    for h = 1:numhist
        IAEMDE = importdata(['FM', int2str(mix), '/n6/IAEHistMDE', int2str(h), '.txt']);
        IAEminTheta(h) = IAEMDE(ind);
    end
    
    minMIAEMDE(3)=mean(IAEminTheta);
    minSDMDE (3)=std(IAEminTheta);
    
    figure
    subplot(2,1,1)
    plot(thetaStar, 'r.')
    hold on
    plot(minIAEMDE, 'b.')
    subplot(2,1,2)
    plot(thetaStar./minIAEMDE, 'b.')
    %============================n=10^6=====================================%
end

%================= n comparisons========================================%
%% optimal bin
figure
color = ['r', 'b', 'g'];

%%PQ vs RegHist
for i = startN:endN
    subplot(2,1,1)
    plot(minBinReg(i), minBinPQ(i), 'o', 'MarkerEdgeColor', color(i), 'MarkerFaceColor', color(i));
    hold on
    title('Number of bins: RegHist vs PQ Hist')
    xlabel('Optimal number of bins for RegHist')
    ylabel('Optimal number of bins for PQ Hist')
    
    subplot(2,1,2)
    plot(minMIAEReg(i), minMIAEPQ(i), 'o', 'MarkerEdgeColor', color(i), 'MarkerFaceColor', color(i));
    hold on
    title('IAE: RegHist vs PQ Hist')
    xlabel('Minimum IAE for RegHist')
    ylabel('Minimum IAE for PQ Hist')
end

subplot(2,1,1)
legend('n=10^4', 'n=10^5', 'n=10^6', 'Location', 'Southeast')
subplot(2,1,2)
legend('n=10^4', 'n=10^5', 'n=10^6', 'Location', 'Southeast')

%%MDE vs PQ
figure
color = ['r', 'b', 'g'];

for i = startN:endN
    subplot(2,1,1)
    plot(minBinMDE(i), minBinPQ(i), 'o', 'MarkerEdgeColor', color(i), 'MarkerFaceColor', color(i));
    hold on
    title('Number of bins: MDE vs PQ Hist')
    xlabel('Optimal number of bins for MDE')
    ylabel('Optimal number of bins for PQ Hist')
    
    subplot(2,1,2)
    plot(minMIAEMDE(i), minMIAEPQ(i), 'o', 'MarkerEdgeColor', color(i), 'MarkerFaceColor', color(i));
    hold on
    title('IAE: MDE vs PQ Hist')
    xlabel('Minimum IAE for MDE')
    ylabel('Minimum IAE for PQ Hist')
end

subplot(2,1,1)
legend('n=10^4', 'n=10^5','n=10^6','Location', 'Southeast')
subplot(2,1,2)
legend('n=10^4', 'n=10^5','n=10^6', 'Location', 'Southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time
figure

% PQ
if (a==1)
    subplot(3,1,1)
    load TimePQ10000
    count = count1/10000;
    timePQ = mean(timePQ);
    plot(count, timePQ, 'r.')
end

hold on

if (b==1)
    subplot(3,1,1)
    load TimePQ100000
    count = 500:500:50000;
    count = count2/100000;
    timePQ = mean(timePQ);
    plot(count, timePQ, 'b.')
end

hold on

if (c==1)
    subplot(3,1,1)
    load TimePQ1000000
    count = 5000:5000:50000;
    count = count3/1000000;
    timePQ = mean(timePQ);
    plot(count, timePQ, 'g.')
end
hold off

title('Time needed to build PQ histogram')
xlabel('Max number of points in each box/n')
ylabel('Time (s)')
legend('n=10^4', 'n=10^5', 'n=10^6')

% MDE
subplot(3,1,2)

if (a==1)
    load TimeMDE10000
    plot(1, timeMDE, 'ro', 'MarkerFaceColor', 'r')
end
hold on
if (b==1)
    load TimeMDE100000
    plot(2, timeMDE, 'bo', 'MarkerFaceColor', 'b')
end
hold on
if (c==1)
    load TimeMDE1000000
    plot(3, timeMDE, 'go', 'MarkerFaceColor', 'g')
end
hold off

title('Time needed to build MDE PQ histogram')
xlabel('n')
ylabel('Time (s)')

% Reg Hist
subplot(3,1,3)
if (a==1)
    load TimeReg10000
    plot(1, timeReg, 'ro', 'MarkerFaceColor', 'r')
end
hold on
if (b==1)
    load TimeReg100000
    plot(2, timeReg, 'bo', 'MarkerFaceColor', 'b')
end
hold on
if (c==1)
    load TimeReg1000000
    plot(3, timeReg, 'go', 'MarkerFaceColor', 'g')
end
hold off

title('Time needed to build regular histogram')
xlabel('n')
ylabel('Time (s)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('==============Output===================')
% container to keep optimal bin size
for i = 1:length(minMIAEPQ)
    i
    [minMIAEPQ(i) minSDPQ(i) minBinPQ(i)]
    [minMIAEReg(i) minSDReg(i) minBinReg(i)]
    [minMIAEMDE(i) minSDMDE(i) minBinMDE(i)]
end


