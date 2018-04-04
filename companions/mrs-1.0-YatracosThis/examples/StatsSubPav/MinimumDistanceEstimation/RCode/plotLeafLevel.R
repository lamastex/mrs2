plotMap1D = function(mydata, str){
  ##Purpose: to plot uniform mixtures given leaf levels
  #Input arguments:
  # mydata: the data (in data frame)
  # str: if str = 1, a plot will be output
  ################################################

function [fhat vol] = plotLeafLevel(leafLevel, emptyLeafLevel, plotting)

if (length(leafLevel) ~= length(emptyLeafLevel) )
  error('Length of 1st arg and length of 2nd arg are different!')
end

%vol
vol = 1./(2.^leafLevel);

Lb = [];
Lb(1) = 0;
for i = 1:(length(vol)-1)
Lb(i+1) = Lb(i) + vol(i);
end
Lb = Lb';

Ub = [];
Ub(1) = vol(1);
for i = 1:(length(vol)-1)
Ub(i+1) = Ub(i) + vol(i+1);
end
Ub = Ub';


%Density estimate
filled = find(emptyLeafLevel == 1);
sizeFilled = length(filled);
fhat = 1/sizeFilled * 2.^leafLevel;
%fhat = 1/length(leafLevel) * 2.^leafLevel;
fhat = fhat';


%empty holes at emptyleaflevels
holes = find(emptyLeafLevel == 0);
fhat(holes) = 0;

if (plotting==1)
YY=[0 1 1 0]' * fhat'; YY=YY(:);
XX=[1 1 0 0]' * Lb' + [0 0 1 1]' *Ub'; XX=XX(:);
plot(XX, YY, 'k', 'LineWidth', 2)
%hold on
%fill(XX, YY, 'c')

elseif (plotting==2)
for i = 1:length(fhat)
plot([Lb(i) Ub(i)], [fhat(i) fhat(i)], 'k', 'LineWidth', 2)
hold on
end

 } #end of plotLeafLevel.R