function [ thetaHat ] = lag1theta(corrVal,dx,corFun)
% Function to estimate theta from lag-1 correlation
%
% Input:
% corrVal - lag-1 correlation
% dx - lag-1 seperation distnace
% corFun - desired correlation function. Options: 'exp','sexp','poly','tri'

    corrVal(corrVal < 0) = 0;

    switch corFun
        case {'exp','markov'}
            thetaHat = -(2*dx)./log(corrVal);
        case {'gauss','sexp'}
            thetaHat = dx*sqrt(pi./-log(corrVal));
        case 'poly'
            thetaHat = dx./(corrVal.^(-1/3)-1);
        case 'tri'
            thetaHat = dx./(1-corrVal);
    end
     
    %if ~isreal(thetaHat)
    %    thetaHat = 0;
    %end


end

