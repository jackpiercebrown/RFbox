function [corrVec,covVec] = corrPairs(samples,mu)
% corrPairs - non-Fourier implementation of sample autocorrelation
%
% Inputs:
% -------
%   samples       - nx by nSamples matrix of data
%   mu (optional) - mean of random field
%
% Outputs:
% --------
%   acf - vector of correlations

% Copyright Jack Pierce-Brown 2018

    [nLags,nSamples] = size(samples);

    if nargin < 2
        mu = mean(mean(samples));
    end

    samplesAdj = samples - mu;

    covVec = zeros(nLags,1);

    divisor = nLags;

    for iLag = 0:(nLags-1)
        sum = 0;
        for iSample = 1:nSamples
            for ix = 1:(nLags - iLag)
                sum = sum + samplesAdj(ix,iSample)*samplesAdj(ix+iLag,iSample);
            end
        end
        covVec(iLag+1) = sum/(nSamples*divisor);
    end
    corrVec = covVec/covVec(1);
    
end

