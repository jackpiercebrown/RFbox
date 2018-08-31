% example.m is a script providing examples of how to use estimation tools.

%% Example 1 - inference from data

% simulate data
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
samples = CMD1D(points,corFun,corLen,dst,mu,sigma,method,nSamples);

% visualise data
fig1 = figure(1);
plot(samples)

% method of moments
mu = mean(samples(:))
sigma = std(samples(:))

% lag-1 correlation
corrVal = lag1corr(samples);
thetaHat = lag1theta(corrVal,dx,corFun)

% acf
ACF2 = corrFFT(samples);
thetaHat = fitACF(ACF2,dx,corFun,lx);

fig2 = figure(2);
plot(ACF2)

% maximum likelihood estimation
[thetaMLE,muMLE,sigmaMLE] = maxLfun(samples,points,corFun,0,lx)

%% Example 2 - Kriging
uncondPoints = linspace(0,8,1e3)';
condPoints = [0.5 1.3 1.5 2.8 4 4.2]';
condVals = [7.7 7.6 7.2 6.3 7.1 6.9]';
mu = 7;
theta = 1.5;
sigma = 0.5;

corFun = 'sexp';

[krige,CIupper,CIlower] = krigeIt(condPoints,condVals,uncondPoints,corFun,mu,sigma,theta);

fig3 = figure(3);
hold on
X = [uncondPoints; flipud(uncondPoints)];
Y = [CIlower; flipud(CIupper)];
h = fill(X',Y',[0.8,0.8,0.8]);
set(h,'EdgeColor','none')
plot(uncondPoints,krige,'color',[1 0 0],'LineWidth',2)
scatter(condPoints,condVals,'MarkerEdgeColor',[1 0 0],'LineWidth',2,'Marker','o','SizeData',50,'LineWidth',2,'MarkerFaceColor',[1,1,1]);
hold off
box on
set(fig3,'pos',[100,100,1200,500])