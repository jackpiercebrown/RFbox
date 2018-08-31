function [condSamples] = CMD1Dcond(uncondPoints,cFun,cLen,dst,mu,sigma,method,nSamples,condPoints,condVals)
% CMD1Dcond - conditional random process simulation
%
% Inputs:
% -------
%   uncondPoints - colum vector of desired simulation points
%   corFun       - desired correlation function:
%                  'markov','gauss','poly','tri','markov2'
%   corLen       - correlation length
%   dst          - marginal distribution: 'normal','lognormal'
%   mu           - mean of random field
%   sigma        - standard deviation of random field
%   method       - decomposition method: 'chol','eig'
%   nSamples     - number of realisations of random field
%   condPoints   - column vector of conditional locations
%   condVals     - column vector of conditional values
%
% Outputs:
% --------
%   condSamples - nPoints by nSamples matrix of samples
%
% Example:
% --------
% samplesMat = CMD1Dcond([0:0.01:9.99]','gauss',2,'lognormal',10,3,'eig',5,[3,5,8]',[15,5,15]');
% plot(samplesMat)

% Copyright Jack Pierce-Brown 2018


nx = length(uncondPoints);

extPoints = [uncondPoints; condPoints];
nCond = length(condPoints);
nUncond = length(uncondPoints);
condInd = [nUncond+1:nUncond+nCond]';

[samples,Cbig] = CMD1D(extPoints,cFun,cLen,dst,mu,sigma,method,nSamples);

if strcmp(dst,'lognormal')
   mu = log(mu)/sqrt(1+(sigma/mu)^2);
   samples = log(samples);
   condVals = log(condVals);
end

% Krige conditional values:
C = Cbig(condInd,condInd);
b = Cbig(condInd,:);
beta = C\b;
condKrige = mu + beta'*(condVals - mu);

condSamples = zeros(nx,nSamples);
for iSample = 1:nSamples
    
    % Krige unconditional values:
    uncondknownPoints = samples(condInd,iSample);
    uncondKrige = mu + beta'*(uncondknownPoints - mu);

    % simulate realisations:
    condSamples(:,iSample) = samples(1:nUncond,iSample) + condKrige(1:nUncond) - uncondKrige(1:nUncond);
end

if strcmp(dst,'lognormal')
    disp('yes')
    condSamples = exp(condSamples);
end

end

