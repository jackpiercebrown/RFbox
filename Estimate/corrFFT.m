function [ acf ] = corrFFT(samplesMat,denomType,mu)
% Adaptation of inbuilt function autocorr to allow for multiple
% realisations, known mean and variable denominator.
%
% Input:
% samplesMat - nx by nSamples matrix of realisations
% denomType (optional) - 'cnst' or 'var'
% mu (optional) - if true mean is known
    
    if nargin < 2
        denomType = 'cnst';
    end
    if nargin < 3
        mu = mean(mean(samplesMat));
    end
    
    [nx,~] = size(samplesMat);

    nFFT = 2^(nextpow2(nx)+1);
    F = fft(samplesMat-mu,nFFT);
    F = F.*conj(F);
    acf = ifft(mean(F,2));
    acf = acf(1:nx);
    if strcmp(denomType,'var')
        acf = (nx./[nx:-1:1]').*acf;
    end

    acf = acf./acf(1);
    %acf = acf./max(abs(acf));
    acf = real(acf);
end