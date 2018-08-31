% FFT1D.m simulates a Gaussian stationary 1D random process with an isotropic
% correlation structure via the fast fourier transform (FFT) method.
%
% Input arguments:
% f - struct detailing field parameters:
%   f.distr - marginal distribution either 'normal' or 'lognormal'
%   mu - mean of marginal distribution
%   sigma - standard deviation of marginal distribution
% 
% c - struct detailing correlation structure:
%   corFun - correlation function. Options: 'exp','gauss','pol.y','tri'
%   corLen - correlation length in x-direction
%
% l - struct detailing length of random field
%   lx - length of field in x-direction
%
% n - struct detailing discretisation of random field
%   nSamples - number of realisations required
%   nx - number of elements in x-direction
%
% Output:
% samplesMat - matrix format: nSamples by nx
%
% a - struct with additional information
%   a.runtime - time to generate realisations
%
% written by Jack Pierce-Brown 2016 (jack.piercebrown@nottingham.ac.uk)

% FFT1D.m simulates a Gaussian stationary 2D random field with an isotropic
% correlation structure via the fast fourier transform (FFT) method.
%
% Written by Jack Pierce-Brown, University of Nottinham, 2016 based on
% simulation and analysis of random fields (Fenton, 1990).

function [ RFreduc ] = FFT1D(lx,nx,corFun,corLen,mu,sigma,nSamples)

% Inputs:
% -------
% lx - length of random process
% nx - number of elements to discretise the random process into
% corFun - target correlation function:
%           options: 'markov','gauss','tri','markov2'
% mu - target mean
% sigma - target standard deviation
% nSamples - desired number of realisations
%
% Outputs:
% --------
% sampleMat



% double field discretisation and length due to periodicity of FFT:
nx2 = nx*2; 
lx2 = lx*2+lx/(nx-1);

%% Determine spectral density discretisation
deltaW = 2*pi*(nx2-1)/(nx2*lx2); % frequency discretisation interval

wx = (0:deltaW:deltaW*(nx2-1))'; % frequency discretisation vector

%% Determine one sided SDF

G0 = sigma^2*corLen/pi;

switch corFun
    case 'markov'
        Gvec = G0./(1+(corLen*wx/2).^2);
    case 'gauss'
        Gvec = G0*exp(-(wx.^2*corLen^2/(4*pi)));
    case 'tri'
        Gvec = G0*(sin(wx*corLen/2)./(wx*corLen/2)).^2;
        Gvec(1) = G0;
    case 'markov2'
        Gvec = 256*G0./((16 + corLen^2*wx.^2).^2);
end

size(Gvec)

%% Determening the standard deviation of coefficients Ak and Bk

Gperiodic = [Gvec; Gvec(1)]; % extend G outside domain

% calculate each term seperately:
G1 = Gperiodic(1:nx+1);
G2 = Gperiodic(nx2+1:-1:nx+1); 
% sum terms according to Fenton (1990):
varAk = 0.25*deltaW*(G1+G2);

% double variance of Ak at 1st and midpoint
varAk([1,nx+1]) = 2*varAk([1 nx+1]);

stdAk = sqrt(varAk); % standard deviation of Ak

stdBk = stdAk;
stdBk([1,nx+1]) = 0; % zero variance of Bk in corners

%% Generate random fields
fieldRealisations = zeros(nx2,nSamples); % initialise realisation matrix

for iSample = 1:nSamples
    Ak = zeros(nx2,1); % initialise Ak, Bk
    Bk = zeros(nx2,1);
    
    Ak(1:nx+1) = randn(nx+1,1).*stdAk(1:nx+1); % generate 1st quadrant
    Bk(1:nx+1) = randn(nx+1,1).*stdBk(1:nx+1);
    
    % symmetry line
    Ak(nx2:-1:nx+2) = Ak(2:nx);
    Bk(nx2:-1:nx+2) = Bk(2:nx); % should be negative?

    realField = real(fft2(Ak - 1i*Bk));

    fieldRealisations(:,iSample)=mu + realField;
end

RFreduc = fieldRealisations(1:nx,:);
