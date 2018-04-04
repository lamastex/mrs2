clc
clear all

d=1;

%est function
A = importdata('EstFunctionAfterNormalized.txt');
A = A.data;

heights = A(:,2);
m = max(heights)
find(heights==max(m))

%true function
mu = zeros(d,1);
sigma = ones(d,1);
p = mvnpdf(mu,mu,sigma');

%difference
diff = p - m 

%get the weights
w = importdata('Weights.txt');
max(w)
find(w==max(w))

%get the data
theData = importdata('Data.txt');


%get the final histogram
B = importdata('FinalHist0.txt');
B = B.data;
hb = B(:,2);
max(hb)
