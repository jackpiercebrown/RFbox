function [corMat] = corMat(points,cFun,cLenX,cLenY)
% corMat - return correlation matrix
%
% Inputs:
% -------
%   points - column vector of spatial locations
%   cFun   - desired correlation function:
%            'markov','gauss','poly','tri','markov2'
%   cLenX  - correlation length in the x-direction
%   cLenY  - correlation length in the y-direction
%
% Outputs:
% --------
%   corMat - nPoints by nPoints correlation matrix

% Copyright Jack Pierce-Brown 2018

matSize = size(points);

if matSize(2) == 1
    distVec = pdist(points);
    adjDist = distVec/cLenX;
elseif matSize(2) == 2
    distX = pdist(points(:,1)); 
    distY = pdist(points(:,2));
    adjDist = sqrt((distX/cLenX).^2 + (distY/(cLenY)).^2);
else
    error('points is invalid size')
end

cVec = corFun(adjDist,cFun);

corMat = squareform(cVec);
[m,~] = size(corMat);
corMat(1:m+1:end) = 1;

end

