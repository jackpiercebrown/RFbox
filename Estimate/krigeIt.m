function [krige,CIupper,CIlower] = krigeIt(condPoints,condVals,uncondPoints,corFun,mu,sigma,theta,CIalpha)
% krigeIt - returns Kriging interpolation and confidence intervals
%
% Inputs:
% -------
% condPoints   - spatial or time location of data
% condVals     - value of process at those points
% uncondPoints - points where prediction is desired
% corFun       - correlation function: 'markov','gauss','tri','poly'
% mu           - mean of random field
% sigma        - standard deviation of random field
% theta        - correlation length of random field
% CIalpha      - (optional) confidence level for confidence intervals
%
% Outputs:
% --------
% krige   - kriging interpolation of data
% CIupper - upper confidence interval
% CIlower - lower confidence interval

% Copyright Jack Pierce-Brown 2018

if nargin < 8
    CIalpha = 0.05;
end
    
pointsExt = [uncondPoints; condPoints];

nCond = length(condPoints);
nUncond = length(uncondPoints);
condInd = [nUncond+1:nUncond+nCond];

corrMat = corMat(pointsExt,corFun,theta);
covMat = sigma^2*corrMat;
C = covMat(condInd,condInd);
b = covMat(condInd,:);
beta = C\b;
krigeExt = mu + beta'*(condVals - mu);

krige = krigeExt(1:nUncond);

errorVarExt = sigma^2 - sum(b.*beta)';
errorVar = errorVarExt(1:nUncond);

CIupper = norminv(CIalpha/2,krige,sqrt(errorVar));
CIlower = norminv(1-CIalpha/2,krige,sqrt(errorVar));

end

