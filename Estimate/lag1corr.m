function [corrVal] = lag1corr(samplesMat,mu)
% lag1corr computes an estimation of the correlation length considering
% only lag-1 seperations.
%
% Input:
% samplesMat - n by nSamples matrix
% mu - (optional argument) if mean of the process is known
%
% Output:
% thetaHat - estimate of correlation length. Note that if computed
% correlation is negative this value is set to zero
% corVal - correlation value from which thetaHat is estimated

    if nargin < 2
        mu = mean(mean(samplesMat));
    end

    [nx,nSamples] = size(samplesMat);

    index = 1;
    nPairs = (nx-1)*nSamples;
    vec1 = zeros(nPairs,1);
    vec2 = zeros(nPairs,1);
    for iSample = 1:nSamples
        for ix = 1:(nx-1)
            vec1(index) = samplesMat(ix,iSample);
            vec2(index) = samplesMat(ix+1,iSample);
            index = index + 1;
        end
    end
    
    vec1adj = vec1 - mu;
    vec2adj = vec2 - mu;
    corrVal = sum(vec1adj.*vec2adj)./((sum(vec1adj.^2)^0.5)*(sum(vec2adj.^2)^0.5));
    
end