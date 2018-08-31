function [ thetaHat ] = lag1theta(corrVal,dx,corFun)
% lag1theta - return correlation length estimate from lag-1 correlation
%
% Inputs:
% -------
% corrVal - lag-1 correlation
% dx      - lag-1 seperation (lag)
% corFun  - desired correlation function: 'markov','gauss','tri','poly'
%
% Outputs:
% --------
% thetaHat - estimate of the correlation length

% Copyright Jack Pierce-Brown 2018

corrVal(corrVal < 0) = 0;

switch corFun
    case {'markov','exp'}
        thetaHat = -(2*dx)./log(corrVal);
    case {'gauss','sexp'}
        thetaHat = dx*sqrt(pi./-log(corrVal));
    case {'poly','polynomial'}
        thetaHat = dx./(corrVal.^(-1/3)-1);
    case {'tri','triangular'}
        thetaHat = dx./(1-corrVal);
end

if ~isreal(thetaHat)
    thetaHat = 0;
end


end

