function [corrVec,covVec] = corrPairs(samples,flag,meanVal)
% corrPairs computes the autocorrelation usiong non-FFT method by
% constructing two vectors and computing sample correlation.
%
% Inputs:
% samples - nx by nSamples matrix of realisations
% meanVal - mean value
% flag - flag for denominator. Default is constant denominator producing
% the same (damped) output as autocorr. If flag == 'varDenom' then ACF is
% more accurate yet unstable for large lags.

    [nLags,nSamples] = size(samples);

    if nargin < 2
        flag = 0;
    end

    if nargin < 3
        meanVal = mean(mean(samples));
    end

    samplesAdj = samples - meanVal;

    covVec = zeros(nLags,1);

    divisor = nLags;

    for iLag = 0:(nLags-1)
        if strcmp(flag,'varDenom')
            divisor = nLags - iLag;
        end
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

