function [fhat, Lb, Ub, vol, count]=plotEstFunction1D(data, what, col, lw)

%function [fhat, Lb, Ub]=subpavhist(data, what, col, lw)
%Calling Function:
%           [fhat, Lb, Ub]=subpavhist(data, what, col, lw)
%
%Purpose: To plot the histogram using subpaving boxes.
%
%Input argument list:
%       data: histogram data
%       what: 1 - need to normalize height; 2 - use heights
%       col: colour
%       lw: line width
%Output argument list:
%       fhat: density estimate
%       Lb  : lower bound
%       Ub  : upper bound
%       vol  : the Lebesgue measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ub and Lb of boxes (output of subpaving should already be in order)
Lb=data(:,3);
Ub=data(:,4);

%Number of points in each box
%this will be the 3rd column of the input

%bin width
vol=data(:,1);

if (what == 1)
    total = sum(data(:,2));
fhat=data(:,2)/total./vol;

else 
    fhat = data(:,2);
end

LL = length(Lb);
newLb = [Lb; Ub(LL)];
newfhat = [fhat; fhat(LL)];
stairs(newLb, newfhat, 'Color', col, 'LineWidth', lw);
    
    YY=[0 1 1 0]' * fhat'; 
    YY=YY(:);
    XX=[1 1 0 0]' * Lb' + [0 0 1 1]' *Ub'; 
    XX=XX(:);
    %plot(XX, YY, 'Color', 'k');
    %fill(XX, YY, [0.5, 0.5, 0.5])
    %fill(XX, YY, 'w')
  
    
    %grid on
    xlabel('support'), ylabel('density')
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of function