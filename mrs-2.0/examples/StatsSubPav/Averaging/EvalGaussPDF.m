function PDF = EvalGaussPDF(x,Mus,CovMat);
% EvalGaussPDF Multivariate normal probability density function.
%
%	PDF = EvalGaussPDF(X,MU,COVM) Returns the value of the multivariate
%	probability density function at the domain locations given in X.
%
%	INPUTS:		X is an n x d matrix of locations in domain
%				MUS is a 1 x d vector
%				COVMAT is a d x d covariance matrix
%
[n,d]=size(x);
PDF = zeros(n,1);
NormConst=(2*pi)^(d/2)*sqrt(det(CovMat));
InvCovMat = inv(CovMat);
for i = 1:n
	Centeredx = x(i,:)-Mus;
	CxICMCxT=Centeredx*InvCovMat*Centeredx';
	PDF(i)=exp((-.5)*CxICMCxT);
end
PDF=PDF/NormConst;
