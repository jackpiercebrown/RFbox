% FFT2D.m simulates a Gaussian stationary 2D random field with an isotropic
% correlation structure via the fast fourier transform (FFT) method.
%
% Written by Jack Pierce-Brown, University of Nottinham, 2016 based on
% simulation and analysis of random fields (Fenton, 1990).

function [ RFreduc ] = FFT2D(lx,ly,nx0,ny0,mu,sigma,corFun,thetaX,thetaY,nSamples)

nx = 2^nextpow2(nx0); % discretisation in x-direction (ensure power of 2)
ny = 2^nextpow2(ny0); % discretisation in y-direction (ensure power of 2)

% double field discretisation and length due to periodicity of FFT:
nx2 = nx*2; 
ny2 = ny*2;
lx2 = lx*2+lx/(nx-1);
ly2 = ly*2+ly/(ny-1);


%% Determine spectral density discretisation
deltaWx = 2*pi*(nx2-1)/(nx2*lx2); % frequenthetaY discretisation interval
deltaWy = 2*pi*(ny2-1)/(ny2*ly2);
deltaW = deltaWx*deltaWy;

wx = 0:deltaWx:deltaWx*(nx2-1); % frequenthetaY discretisation vector
wy = 0:deltaWy:deltaWy*(ny2-1);

%% Determine one sided SDF
[wxMat,wyMat] = meshgrid(wx,wy);
if strcmp(corFun,'markov')==1 % exponential (Markov) spectral density function
    Gmat = (sigma^2*thetaX*thetaY)./(2*pi*(1+(thetaX*wxMat./2).^2 + (thetaY*wyMat./2).^2).^(3/2))'; 
elseif strcmp(corFun,'gauss')==1 % square-exponential (Gaussian) sdf
    Gmat = (sigma^2*thetaX*thetaY/pi^2)*exp(-(wxMat.^2*thetaX^2+wyMat.^2*thetaY^2)./(4*pi))';
end

%% Determening the standard deviation of coefficients Ak and Bk

Gperiodic = [Gmat Gmat(:,1); Gmat(1,:) Gmat(1,1)]; % extend G outside domain

% calculate each term seperately:
G1 = Gperiodic(1:nx+1,1:ny+1);
G2 = Gperiodic(nx2+1:-1:nx+1,1:ny+1);
G3 = Gperiodic(1:nx+1,ny2+1:-1:ny+1);
G4 = Gperiodic(nx2+1:-1:nx+1,ny2+1:-1:ny+1);

% sum terms according to Fenton (1990):
varAk = 0.125*deltaW*(G1+G2+G3+G4);

% double variance of Ak in quadrant corners:
varAk([1 nx+1],[1 ny+1]) = 2*varAk([1 nx+1],[1 ny+1]);

stdAk = sqrt(varAk); % standard deviation of Ak

stdBk = stdAk;
stdBk([1 nx+1],[1 ny+1]) = 0; % zero variance of Bk in corners

%% Generate random fields
fieldRealisations = zeros(nx2,ny2,nSamples); % initialise realisation matrix

for iSample = 1:nSamples
    Ak = zeros(nx2,ny2); % initialise Ak, Bk
    Bk = zeros(nx2,ny2);
    
    Ak(1:nx+1,1:ny+1) = randn(nx+1,ny+1).*stdAk(1:nx+1,:); % generate 1st quadrant
    Bk(1:nx+1,1:ny+1) = randn(nx+1,ny+1).*stdBk(1:nx+1,:);
    
    % symmetry line 1
    Ak(nx2:-1:nx+2,1:ny+1) = Ak(2:nx,1:ny+1);
    Bk(nx2:-1:nx+2,1:ny+1) = Bk(2:nx,1:ny+1); % should be negative?
    
    % symmetry line 2
    Ak(1:nx2,ny2:-1:ny+2) = Ak(1:nx2,2:ny);
    Bk(1:nx2,ny2:-1:ny+2) = -Bk(1:nx2,2:ny);
   
    realField = real(fft2(Ak - 1i*Bk));

    fieldRealisations(:,:,iSample)=mu + realField;
end

RFreduc = fieldRealisations(1:nx,1:ny)';