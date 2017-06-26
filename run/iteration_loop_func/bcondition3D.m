function BBV = bcondition3D(BBvector)

% In this function the dirichlet boundary condition is applied at the
% control points located at the boundary of the image domain
% The end condition is fixed, i.e. the displacement of the control points
% at the boundary is zero.

% INPUT:
% BBvector: the RHS matrix containing the the value of displacement of the
% control points, x,y,z directions. The fourth index corresponds to whether
% the control point lies on the boundary of the image domain (1) or not
% (0).

% OUTPUT:
% BBV: the output RHS matrix after applying the Dirichlet BC

%%

index = find(BBvector(:,4)==1);

BBvector1 = BBvector;

vec1 = zeros(size(index,1),3);
BBvector1(index,1:3) = vec1;


BBV = BBvector1;

end