function [thetaHat,fittedACF] = fitACF(ACF,dx,corFun,theta0,weights,tol)
% fitACF - estimate correlation length from sample correlation function
%
% Inputs:
% -------
%   ACF     - sample correlation function
%   dx      - spatial seperation (lag size)
%   cFun    - desired correlation function:
%             'markov','gauss','poly','tri','markov2'
%   theta0  - start point of optimisation algorithm
%   weights - weights for regression
%   tol     - stopping tolerance for convergence
%
% Ouputs:
% -------
% thetaHat  - estimated value of correlation length
% fittedACF - fit theoretical correlation function to data

% Copyright Jack Pierce-Brown 2018

[nLag,nACF] = size(ACF);

if nargin < 5
    weights = ones(nLag,1);
end
if nargin < 6
    tol = 1e-12;
end

nW = length(weights);
if nW < nLag
    weightPad = zeros(nLag,1);
    weightPad(1:nW) = weights;
    weights = weightPad;
end

lagVec = [0:dx:(nLag-1)*dx]';
lagMat = repmat(lagVec,1,nACF);

switch corFun
    case {'markov','exp'}
        f = @(theta) exp(-2*abs(lagMat)./theta);
    case {'gauss','sexp'}
        f = @(theta) exp(-pi*(abs(lagMat)./theta).^2);
    case {'tri','triangular'}
        f = @(theta) (1 - lagMat/theta).*[ones(theta/dx,1);zeros(nLag-theta/dx,1)];
end

weightMat = repmat(weights,[1,nACF]);

errorF = @(theta) weightMat.*(f(theta) - ACF).^2;
thetaHat = lsqnonlin(errorF,theta0,[],[],optimset('Display','off','TolFun',tol));

fittedACF = f(thetaHat);
    
end

