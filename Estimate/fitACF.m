function [thetaHat,fittedACF] = fitACF(acf,dx,corFun,theta0,weights,tol)
% fitACF fits a theoretical correlation function to autocorrelation data
% using non-linear least-squares LSQNONLIN to minimise squared error
%
% Inputs:
% acf - autocorrelation data (single ACF in column vector or many in
% matrix)
% corFun - theoretical correlation model to be fitted
% weights - weighting function
%
% Outputs:
% thetaHat - value of theta that minimises misfit

    [nLag,nACF] = size(acf);

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
        case 'sexp'
            f = @(theta) exp(-pi*(abs(lagMat)./theta).^2);
        case 'tri'
            f = @(theta) (1 - lagMat/theta).*[ones(theta/dx,1);zeros(nLag-theta/dx,1)];
    end
    
    weightMat = repmat(weights,[1,nACF]);
  
    errorF = @(theta) weightMat.*(f(theta) - acf).^2;
    thetaHat = lsqnonlin(errorF,theta0,[],[],optimset('Display','off','TolFun',tol));
    
    fittedACF = f(thetaHat);
    
end

