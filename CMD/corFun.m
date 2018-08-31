function [corVec] = corFun(distVec,corFun)
% corFun - returns correlation vector for desired correlation function
%
% Inputs:
% -------
%   distVec - vector of spatial or time seperations
%   corFun  - desired correlation function:
%            'markov','gauss','poly','tri','markov2'
%
% Outputs:
% --------
%   corVec - column vector of correlations

% Copyright Jack Pierce-Brown 2018

switch corFun
    case {'markov','exp'}
        corVec = exp(-2*distVec);
    case {'gauss','sexp'}
        corVec = exp(-pi*distVec.^2);
    case {'poly','polynomial'}
        corVec = (1 + distVec).^-3;
    case {'tri','triangular'}
        corVec = 1 - distVec;
        corVec(corVec < 0) = 0;
    case {'markov2','smarkov'};
        corVec = (1 + 4*distVec).*exp(-4*distVec);
    otherwise
        error('Invalid correlation function')
end

end

