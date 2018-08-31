function [samples] = FFT1D(lx,nx0,cFun,cLen,mu,sigma,nSamples)
% FFT1D - simulate random process using fast Fourier transform method
%
% Inputs:
% -------
%   lx       - length of random process
%   nx0      - number of elements in x-direction
%   cFun     - correlation function:
%              'markov','gauss','tri','markov2'
%   mu       - mean of random field
%   sigma    - standard deviation of random process
%   nSamples - number of realisations of random process
%
% Outputs:
% --------
%   sampleMat - nx by nSamples matrix of realisations
%
% Example:
% --------
% samples = FFT1D(10,1e3,'markov',2,10,3,5);
% plot(samples)

% Copyright Jack Pierce-Brown 2018

% double field due to periodicity of FFT:
nx = 2^nextpow2(nx0);
nx2 = nx*2; 
lx2 = lx*2+lx/(nx-1);

% discretise frequency spectrum:
deltaW = 2*pi*(nx2-1)/(nx2*lx2);
wx = (0:deltaW:deltaW*(nx2-1))';

G0 = sigma^2*cLen/pi;

switch cFun
    case 'markov'
        Gvec = G0./(1+(cLen*wx/2).^2);
    case 'gauss'
        Gvec = G0*exp(-(wx.^2*cLen^2/(4*pi)));
    case 'tri'
        Gvec = G0*(sin(wx*cLen/2)./(wx*cLen/2)).^2;
        Gvec(1) = G0;
    case 'markov2'
        Gvec = 256*G0./((16 + cLen^2*wx.^2).^2);
end

% determine standard deviation of fourier coefficients Ak and Bk
Gperiodic = [Gvec; Gvec(1)];
G1 = Gperiodic(1:nx+1);
G2 = Gperiodic(nx2+1:-1:nx+1); 
varAk = 0.25*deltaW*(G1+G2);
varAk(nx+1) = 2*varAk(nx+1); % double variance of Ak at midpoint

stdAk = sqrt(varAk);
stdBk = stdAk;
stdBk([1,nx+1]) = 0; % zero variance of Bk in corners

% simulate realisatrions
fieldRealisations = zeros(nx2,nSamples);

for iSample = 1:nSamples
    Ak = zeros(nx2,1); 
    Bk = zeros(nx2,1);
    
    % generate 1st half
    Ak(1:nx+1) = randn(nx+1,1).*stdAk(1:nx+1);
    Bk(1:nx+1) = randn(nx+1,1).*stdBk(1:nx+1);
    
    % 2nd half by symmetry
    Ak(nx2:-1:nx+2) = Ak(2:nx);
    Bk(nx2:-1:nx+2) = -Bk(2:nx);

    realField = real(fft2(Ak - 1i*Bk));
    fieldRealisations(:,iSample) = mu + realField;
end

samples = fieldRealisations(1:nx,:);

end
