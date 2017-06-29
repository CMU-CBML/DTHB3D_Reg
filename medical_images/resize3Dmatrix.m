function [I ] = resize3Dmatrix( nx,ny,nz,I )
%source code taken from the website https://github.com/stellaccl/cdmffd-image-registration
%Paper: Chan, Chiu Ling, et al. "Two and Three Dimensional Image Registration Based on B-Spline Composition
% and Level Sets." Communications in Computational Physics 21.2 (2017): 600-622.
% resize3Dmatrix resize 3 dimensional matrix (I) to a given size (nx,ny,nz)
% Input: I = 3 dimensional matrix. 
%        nx, ny,nz = desired size in x, y ,and z direction.
% Output: I = rezised 3 dimensional matrix. 

[y, x, z]=...
   ndgrid(linspace(1,size(I,1),ny),...
          linspace(1,size(I,2),nx),...
          linspace(1,size(I,3),nz));
I=interp3(I,x,y,z);
end

