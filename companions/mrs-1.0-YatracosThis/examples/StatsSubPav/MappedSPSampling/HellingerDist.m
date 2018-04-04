clc
clear all
close all

%get hellinger distance
d=1;

%import actual data
actual = importdata('ActualData1.txt');
mActual = mean(actual);
vActual = cov(actual);

mMapped = zeros(d,1)';
vMapped = eye(d,d);
 
%import mapped data
% mapped = importdata('MappedData1.txt');
% mMapped = mean(mapped);
% vMapped = cov(mapped);

%first get the difference between the mean vectors
diffMean = mActual - mMapped;


%for 1D
if (d==1)
    HD1 = 1 - sqrt(2*sqrt(vMapped)*sqrt(vActual)/(vMapped+vActual))* ...
                exp(-0.25*(diffMean)^2/(vMapped+vActual));
    HD1 = sqrt(HD1)
    
else
    %get matrix P where P = (vActual+vMapped)/2
    P = (vActual+vMapped)/2;
    inv(P)
    
    (diffMean) * inv(P) * (diffMean)'
    
    D = 1/8 *(diffMean) * inv(P) * (diffMean)' + ...
        0.5*log(det(P)/sqrt(det(vActual) * det(vMapped)));
    
    det(vMapped)
    det(vActual)
    det(P)
    
   
    
    BC = exp(-1*D);
    
   
    
    HD = sqrt(1-BC)
end