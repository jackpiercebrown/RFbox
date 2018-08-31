function [samples] = FFT2D(lx,ly,nx0,ny0,mu,sigma,corFun,thetaX,thetaY,nSamples)
% FFT2D - simulate random process using fast Fourier transform method
%
% Inputs:
% -------
%   lx       - length of random process
%   nx0      - number of elements in x-direction
%   ny0      - number of elements in y-direction
%   cFun     - correlation function:
%              'markov','gauss','tri','markov2'
%   mu       - mean of random field
%   sigma    - standard deviation of random field
%   nSamples - number of realisations of random field
%
% Outputs:
% --------
%   sampleMat - nx by nSamples matrix of realisations
%
% Example:
% --------
% samples = FFT2D(10,8,1e3,800,10,3,'markov',2,3,5);
% image(samples(:,:,1),'CDataMapping','scaled')

% Copyright Jack Pierce-Brown 2018

% ensure power of 2 and double field size in each dimension:
nx = 2^nextpow2(nx0);
ny = 2^nextpow2(ny0);
nx2 = nx*2; 
ny2 = ny*2;
lx2 = lx*2+lx/(nx-1);
ly2 = ly*2+ly/(ny-1);

% discretise frequency spectrum:
deltaWx = 2*pi*(nx2-1)/(nx2*lx2);
deltaWy = 2*pi*(ny2-1)/(ny2*ly2);
deltaW = deltaWx*deltaWy;
wx = 0:deltaWx:deltaWx*(nx2-1);
wy = 0:deltaWy:deltaWy*(ny2-1);

[wxMat,wyMat] = meshgrid(wx,wy);
if strcmp(corFun,'markov') == 1
    Gmat = (sigma^2*thetaX*thetaY)./(2*pi*(1+(thetaX*wxMat./2).^2 + (thetaY*wyMat./2).^2).^(3/2))'; 
elseif strcmp(corFun,'gauss') == 1
    Gmat = (sigma^2*thetaX*thetaY/pi^2)*exp(-(wxMat.^2*thetaX^2+wyMat.^2*thetaY^2)./(4*pi))';
end

% determine standard deviation of fourier coefficients Ak and Bk
Gperiodic = [Gmat Gmat(:,1); Gmat(1,:) Gmat(1,1)];
G1 = Gperiodic(1:nx+1,1:ny+1);
G2 = Gperiodic(nx2+1:-1:nx+1,1:ny+1);
G3 = Gperiodic(1:nx+1,ny2+1:-1:ny+1);
G4 = Gperiodic(nx2+1:-1:nx+1,ny2+1:-1:ny+1);
varAk = 0.125*deltaW*(G1+G2+G3+G4);

% double variance of Ak in quadrant corners:
varAk([1 nx+1],[1 ny+1]) = 2*varAk([1 nx+1],[1 ny+1]);
stdAk = sqrt(varAk);
stdBk = stdAk;
stdBk([1 nx+1],[1 ny+1]) = 0; % zero variance of Bk in corners

% simulate realisatrions
fieldRealisations = zeros(nx2,ny2,nSamples);

for iSample = 1:nSamples
    Ak = zeros(nx2,ny2);
    Bk = zeros(nx2,ny2);
    
    % generate first quadarant
    Ak(1:nx+1,1:ny+1) = randn(nx+1,ny+1).*stdAk(1:nx+1,:); 
    Bk(1:nx+1,1:ny+1) = randn(nx+1,ny+1).*stdBk(1:nx+1,:);
    
    % symmetry line 1
    Ak(nx2:-1:nx+2,1:ny+1) = Ak(2:nx,1:ny+1);
    Bk(nx2:-1:nx+2,1:ny+1) = Bk(2:nx,1:ny+1); % should be negative?
    
    % symmetry line 2
    Ak(1:nx2,ny2:-1:ny+2) = Ak(1:nx2,2:ny);
    Bk(1:nx2,ny2:-1:ny+2) = -Bk(1:nx2,2:ny);
   
    realField = real(fft2(Ak - 1i*Bk));
    fieldRealisations(:,:,iSample) = mu + realField;
end

samples = fieldRealisations(1:nx0,1:ny0,:);

end