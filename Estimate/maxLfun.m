function [thetaMLE,muMLE,sigmaMLE,lval] = maxLfun(samples,points,corFun,lowerTheta,upperTheta)
% maxLfun - returns parameter estimates that maximise Gaussian likelihood
%
% Inputs:
% -------
% samples    - nx by nSamples matrix of data
% points     - spatial locations of samples
% corFun     - correlation function: 'markov','gauss','tri','poly'
% lowerTheta - lower bound on the correlation length
% upperTheta - upper bound on the correlation length
%
% Outputs:
% --------
% thetaMLE - estimate of the correlation length
% muMLE    - estimate of the mean
% sigmaMLE - estimate of the standard deviation
% lval     - value of maximum likelihood

% Copyright Jack Pierce-Brown 2018

f = @(theta) -gaussLfun(samples,corMat(points,corFun,theta));

thetaMLE = fminbnd(f,lowerTheta,upperTheta);
[lval,muMLE,sigmaMLE] = gaussLfun(samples,corMat(points,corFun,thetaMLE));

if thetaMLE == upperTheta
    printf('Warning: consider increasing upperbound')
end

end