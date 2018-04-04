%% get MIAE and time

totalIAE = [];
totalTime =[];

for i = 1:23

    %filename = ['UnifIAE', int2str(i), 'txt'];
    %filename = ['GaussianIAE', int2str(i), 'txt'];
    filename = ['GaussianIAE', int2str(i-1), '.txt'];
    
    %timeF = ['UnifTime', int2str(i), 'txt'];
    %timeF = ['GaussianTime', int2str(i), 'txt'];
    timeF = ['GaussianTime', int2str(i-1), '.txt'];
    
    IAE = importdata(filename);
    totalIAE(i) = IAE;
    
    times = importdata(timeF);
    totalTime(i) = times;
end

mean(totalIAE)
std(totalIAE)

mean(totalTime)
std(totalTime)




