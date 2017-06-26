function [ct1, ct2] = computeCoeffMat(pU, pV, pW,ActiveFlagF,NewInactiveF,cmat, local_s1)
%#codegen
% This function finds the coefficients of the active splines which are non
% zero over the active element

% INPUT:
% pU: degree of spline in u direction
% pV: degree of spline in v direction
% pW: degree of spline in w direction
% ActiveFlagF: the active flag matrix of next refinement level
% NewInactiveF: the newly de-activated flag matrix of next level
% cmat: coefficient matrix assembled 
% local_s1: non zero splines over the element

% OUTPUT:
% ct1: old coefficient matrix
% ct2: new coefficient matrix
ct = zeros((pU+1)*(pV+1)*(pW+1));
ct11 = zeros((pU+1)*(pV+1)*(pW+1));

for j = 1:(pU+1)*(pV+1)*(pW+1)
    [~, Icmat, Ilocal] = intersect(cmat(j,:,2),local_s1,'stable');
    ct(j,Ilocal) = cmat(j,Icmat,1);

    indx1 = find((ActiveFlagF(:,1)==0 & NewInactiveF(:,1)==0));
    ct11(j,indx1) = ct(j,indx1);
    
end

ct1 = ct;
ct2 = ct11;

end