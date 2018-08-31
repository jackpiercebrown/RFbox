function [thetaMLE,muMLE,sigmaMLE,lval] = maxLfun(samples,points,corFun,lowerTheta,upperTheta)
% This function finds the value of theta that maximises the likelihood 
% function and the associated MLE values of mu and sigma.
%
% Input:
% samples - nx by nSamples matrix
% points - spatial locations of samples
% corFun - desired correlation function: 'exp','sexp','tri'
%

    f = @(theta) -gaussLfun(samples,calcCorrMat(points,corFun,theta));

    thetaMLE = fminbnd(f,lowerTheta,upperTheta);
    [lval,muMLE,sigmaMLE] = gaussLfun(samples,calcCorrMat(points,corFun,thetaMLE));
    
    if thetaMLE == upperTheta
        printf('Consider increasing upperbound')
    end
    
end