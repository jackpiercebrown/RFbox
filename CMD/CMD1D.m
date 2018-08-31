function [samples,cMat] = CMD1D(points,cFun,cLen,dst,mu,sigma,method,nSamples)
% CMD1D - random process simulation through covariance matrix decomposition
%
% Inputs:
% -------
%   points   - column vector of spatial locations
%   cFun     - desired correlation function:
%              'markov','gauss','poly','tri','markov2'
%   cLen     - correlation length
%   dst      - marginal distribution: 'normal','lognormal'
%   mu       - mean of random field
%   sigma    - standard deviation of random field
%   method   - decomposition method: 'chol','eig'
%   nSamples - number of realisations of random field
%
% Outputs:
% --------
%   samplesMat - number of points by nSamples matrix of realisations
%   corrMat    - correlation matrix
%
% Example:
% --------
% samplesMat = CMD1D([0:0.01:9.99]','markov',2,'lognormal',10,1,'chol',5);
% plot(samplesMat)

% Copyright Jack Pierce-Brown 2018

cMat = corMat(points,cFun,cLen);

nx = length(cMat);

if strcmp(dst,'lognormal')
    muNorm = log(mu/sqrt(1+(sigma/mu)^2));
    sigmaNorm = sqrt(log(1 + (sigma/mu)^2));
    cMat = log((exp(sigmaNorm^2)-1)*cMat+1)./sigmaNorm^2;
end

if strcmp(method,'chol')
    decompCorrMat = chol(cMat,'upper');
elseif strcmp(method,'eig')
    [eigVec,eigVal] = eig(cMat);
    if any(eigVal(:) < 0)
        disp('Warning: covariance matrix is not positive semidefinite')
    end
    decompCorrMat = sqrt(abs(eigVal))*eigVec';
end
    
samplesMat = normrnd(0,1,[nSamples,nx])*decompCorrMat;

if strcmp(dst,'normal')
    samples = samplesMat'*sigma + mu;
elseif strcmp(dst,'lognormal')
    samples = exp(samplesMat'*sigmaNorm + muNorm);
end

end

