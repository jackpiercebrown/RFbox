% example.m is a script providing examples of how to use estimation tools.

nx = 1e2;
lx = 10;
dx  = lx/nx;
points = linspace(0,10,nx)';
corFun = 'markov';
corLen = 2;
dst = 'normal';
mu  = 10;
sigma = 3;
method = 'eig';
nSamples = 6;

samples = CMD1D(points,corFun,corLen,dist,mu,sigma,method,nSamples);

%condPoints = [0.5 1.3 1.5 2.8 4 4.2]';
%condVals = [7.7 7.6 7.2 6.3 7.1 6.9]';

% visualise data
fig1 = figure(1);
plot(samples)

%% Method of moments

% mu
mu = mean(samples(:))

% sigma
sigma = std(samples(:))

% lag1
corrVal = lag1corr(samples);
thetaHat = lag1theta(corrVal,dx,corFun)

% full acf (FFT and non-FFT)
%ACF = autocorr(samplesMat(:,1),nx-1);
ACF2 = corrFFT(samples);
figure(2)
plot(ACF2)

thetaHat = fitACF(ACF2,dx,corFun,lx)

% < figure of ACF >

%% Maximum likelihood
% mu, sigma theta
[thetaMLE,muMLE,sigmaMLE] = maxLfun(samples,points,corFun,0,lx);

%% Kriging

condPoints = [1,2,3,4,5,6,7]';
condVals = [5,6,7,6,4,5,4]';

corFun = 'sexp';

[krige,CIupper,CIlower] = krigeIt(condPoints,condVals,uncondPoints,corFun,mu,sigma,theta);

figure(3)
hold on
X = [uncondPoints; flipud(uncondPoints)];
Y = [CIlower; flipud(CIupper)];
h = fill(X',Y',[0.8,0.8,0.8]);
set(h,'EdgeColor','none')
plot(uncondPoints,krige,'color',[1 0 0],'LineWidth',2)
scatter(condPoints,condVals,'MarkerEdgeColor',[1 0 0],'LineWidth',2,'Marker','o','SizeData',50,'LineWidth',2,'MarkerFaceColor',[1,1,1]);
hold off
grid on
box on