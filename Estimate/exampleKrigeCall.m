% script to plot krige and confidence intervals for data

condPoints = [1,2,3,4,5,6,7]';
condVals = [5,6,7,6,4,5,4]';

lx = 10;
nx = 1e3;
dx = lx/nx;
uncondPoints = [dx/2:dx:lx]';

corFun = 'sexp';
lowerTheta = 0;
upperTheta = 100;

[theta,mu,sigma,lval] = maxLfun(condVals,condPoints,corFun,lowerTheta,upperTheta);

[krige,CIupper,CIlower] = krigeIt(condPoints,condVals,uncondPoints,corFun,mu,sigma,theta);

figure(1)
scatter(condPoints,condVals)
hold on
plot(uncondPoints,krige)


fig2 = figure(2);
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
