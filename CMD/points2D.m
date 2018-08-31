function [points] = points2D(lx,ly,nx,ny)
% points2D - returns coordinates of points in 2D mesh
%
% Inputs:
% -------
%   lx - length in x-direction
%   ly - length in y-direction
%   nx - number of elements in x-direction
%   ny - number of elements in y-direction
%
% Outputs:
%   points - 2 by nPoints column vector of spatial coordinates

% Copyright Jack Pierce-Brown 2018

dx = lx/nx;
dy = ly/ny;

points = zeros(nx*ny,2);
iPoint = 1;
for ix = 1:nx
    xVal = ix*dx - dx/2;
    for iy = 1:ny
        yVal = iy*dy - dy/2;
        points(iPoint,:) = [xVal -yVal];
        iPoint = iPoint + 1;
    end
end