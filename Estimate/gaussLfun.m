function [lval,muPartial,sigmaPartial] = gaussLfun(samples,corrMat)
% gaussLfun - returns likelihood for multivariate Gaussian samples
%
% Inputs:
% -------
% samples - nx by nSamples matrix of realisations of a Gaussian process
% corrMat - correlation matrix
%
% Outputs:
% --------
% lval - likelihood value from evaluation of the likelihood function
% muPartial - MLE value of mu for given correlation matrix
% sigmaPartial - MLE value of sigma for given correlation matrix
%
% [1] Fenton (1999) - Estimation for stochastic soil models

% Copyright Jack Pierce-Brown 2018
    
[nx,nSamples] = size(samples);

oneVec = ones(nx,1);
oneVecLarge = repmat(oneVec,nSamples,1);

s = corrMat\oneVec;
sLarge = repmat(s,nSamples,1);
rLarge = zeros(nx*nSamples,1);

for iSample = 1:nSamples 
    r = corrMat\samples(:,iSample);
    ind2 = nx*iSample;
    ind1 = ind2 - (nx - 1);
    rLarge(ind1:ind2) = r;
end

flatSamples = reshape(samples,[numel(samples),1]);

logDet = nSamples*2*sum(log(diag(chol(corrMat))));
muPartial = (oneVecLarge'*rLarge)/(oneVecLarge'*sLarge);
varPartial = (1/(nx*nSamples))*(flatSamples - muPartial*oneVecLarge)'*rLarge;
lval = -(nx*nSamples/2)*log(varPartial) - 0.5*logDet;

sigmaPartial = sqrt(varPartial);

end

