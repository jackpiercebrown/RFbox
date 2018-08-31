% example.m is a script providing examples of how to simulate random fields
% using the fast Fourier transform (FFT) method.

%% Example 1 - 1D random processes

lx = 1;
nx = 1e3;
corLen = 0.2;
mu = 10;
sigma = 3;
nSamples = 1;
fieldMarkov = FFT1D(lx,nx,'markov',corLen,mu,sigma,nSamples);
fieldGauss = FFT1D(lx,nx,'gauss',corLen,mu,sigma,nSamples);
fieldTri = FFT1D(lx,nx,'tri',corLen,mu,sigma,nSamples);
fieldMarkov2 = FFT1D(lx,nx,'markov2',corLen,mu,sigma,nSamples);

fig1 = figure(1);
plot(fieldMarkov)
hold on
plot(fieldGauss)
plot(fieldTri)
plot(fieldMarkov2)
hold off
set(fig1,'pos',[100,100,1200,500])

%% Example 2 - 2D random fields

lx = 1;
ly = 1;
nx0 = 1e3;
ny0 = 800;
mu = 10;
sigma = 0.3;
corFun = 'gauss';
thetaX = 0.5*lx;
thetaY = 0.1*ly;
nSamples = 1;
sampleMarkov = FFT2D(lx,ly,nx0,ny0,mu,sigma,'markov',thetaX,thetaY,nSamples);
sampleGauss = FFT2D(lx,ly,nx0,ny0,mu,sigma,'gauss',thetaX,thetaY,nSamples);

fig2 = figure(2);
image([0 lx],[0 ly],sampleMarkov(:,:,1),'CDataMapping','scaled')
set(fig2,'pos',[100,100,1200,500])

fig3 = figure(3);
image([0 lx],[0 ly],sampleGauss(:,:,1),'CDataMapping','scaled')
set(fig3,'pos',[100,100,1200,500])