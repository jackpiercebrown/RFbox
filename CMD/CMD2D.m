function [samples,samples2D,cMat] = CMD2D(lx,ly,nx,ny,corFun,corLenX,corLenY,dist,mu,sigma,method,nSamples)
% CMD2D - Random field simulation through covariance matrix decomposition
%
% Inputs:
% -------
%   lx       - length of field in x-direction
%   ly       - length of field in y-direction
%   nx       - number of elements in x-direction
%   ny       - number of elements in y-direction
%   cFun     - desired correlation function: 'markov','gauss'
%   cLen     - correlation length
%   dst      - marginal distribution: 'normal','lognormal'
%   mu       - mean of random field
%   sigma    - standard deviation of random field
%   method   - decomposition method: 'chol','eig'
%   nSamples - number of realisations of random field
%
% Outputs:
% --------
%   samples   - nx*ny by nSamples matrix of realisations
%   samples2D - nx by ny by nSamples matrix of realisations
%   cMat      - correlation matrix
%
% Example:
% --------
% [~,fields] = CMD2D(10,8,60,40,'markov',2,3,'normal',10,3,'chol',5);
% image(fields(:,:,1),'CDataMapping','scaled')

% Copyright Jack Pierce-Brown 2018

points = points2D(lx,ly,nx,ny);
nEl = nx*ny;

cMat = corMat(points,corFun,corLenX,corLenY);

if strcmp(dist,'lognormal')
    muNorm = log(mu/sqrt(1+(sigma/mu)^2));
    sigmaNorm = sqrt(log(1 + (sigma/mu)^2));
    cMat = log((exp(sigmaNorm^2)-1)*cMat+1)./sigmaNorm^2;
end

if strcmp(method,'chol')
    decompCorrMat = chol(cMat,'upper');
elseif strcmp(method,'eig')
    [eigVec,eigVal] = eig(cMat);
    decompCorrMat = sqrt(abs(eigVal))*eigVec';
    if any(eigVal(:) < 0)
        disp('Warning: covariance matrix is not positive semidefinite')
    end
end
    
samplesMat = normrnd(0,1,[nSamples,nEl])*decompCorrMat;

if strcmp(dist,'normal')
    samples = samplesMat'*sigma + mu;
elseif strcmp(dist,'lognormal')
    samples = exp(samplesMat'*sigmaNorm + muNorm);
end

samples2D = reshape(samples,[ny, nx, nSamples]);

end