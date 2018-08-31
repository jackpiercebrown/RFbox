function [ acf ] = corrFFT(samples,mu)
% corrFFT - adaptation of autocorr for multiple realisations
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
    
if nargin < 2
    mu = mean(mean(samples));
end

[nx,~] = size(samples);

nFFT = 2^(nextpow2(nx)+1);
F = fft(samples-mu,nFFT);
F = F.*conj(F);
acf = ifft(mean(F,2));
acf = acf(1:nx);

acf = acf./acf(1);
acf = real(acf);

end