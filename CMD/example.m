% example.m is a script providing examples of how to generate random fields
% using the covariance matrix decomposition (CMD) method.

%% Example 1 - 1D random process simulation

points = linspace(0,1,256)';
RFsamples = CMD1D(points,'markov',0.2,'normal',1,0.3,'chol',3);

fig1 = figure(1);
plot(points,RFsamples)
set(fig1,'pos',[100,100,1200,500])

%% Example 2 - 1D random process simulation (more correlation functions)

nx = 512;
lx = 10;
points = linspace(0,lx,nx)';
cLen = 2;
mu = 10;
sigma = 2;
method = 'eig';
nSamples = 1;
RFmarkov = CMD1D(points,'markov',cLen,'normal',mu,sigma,method,nSamples);
RFgauss = CMD1D(points,'gauss',cLen,'normal',mu,sigma,method,nSamples);
RFtri = CMD1D(points,'tri',cLen,'lognormal',mu,sigma,method,nSamples);
RFmarkov2 = CMD1D(points,'markov2',cLen,'lognormal',mu,sigma,method,nSamples);

fig2 = figure(2);
plot(points,RFmarkov)
hold on
plot(points,RFgauss)
plot(points,RFtri)
plot(points,RFmarkov2)
hold off
legend('markov','gauss','tri','markov2')
set(fig2,'pos',[100,100,1200,500])

%% Example 3 - Simulation of 2D random fields

lx = 10;
ly = 5;
nx = 100;
ny = 30;
cFun = 'markov';
cLenX = 1;
cLenY = 2;
dst = 'normal';
mu = 10;
sigma = 1;
method = 'chol';
nSamples = 1e3;
[samplesMat,samples2D] = CMD2D(lx,ly,nx,ny,cFun,cLenX,cLenY,dst,mu,sigma,method,nSamples);

fig3 = figure(3);
image(samples2D(:,:,1),'CDataMapping','scaled')
set(fig3,'pos',[100,100,1200,500])

%% Example 4 - Conditional simulation of 1D random processes
condPoints = [0.5 1.3 1.5 2.8 4 4.2]';
condVals = [7.7 7.6 7.2 6.3 7.1 6.9]';
sigma = 0.5;
mu = 7;
lx = 10;
nx = 1e3;
uncondPoints = linspace(0,lx,nx)';
cFun = 'gauss';
theta = 1.5;
dst = 'lognormal';
method = 'eig';
nSamples = 3;
condSamples = CMD1Dcond(uncondPoints,cFun,theta,dst,mu,sigma,method,nSamples,condPoints,condVals);

fig4 = figure(4);
plot(uncondPoints,condSamples)
hold on
scatter(condPoints,condVals,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1])
hold off
set(fig4,'pos',[100,100,1200,500])