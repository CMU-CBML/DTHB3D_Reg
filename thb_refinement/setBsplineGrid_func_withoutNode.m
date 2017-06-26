function [Basis,Elem] = setBsplineGrid_func_withoutNode(knotvectorU, knotvectorV, knotvectorW,param)
%This function stores the element and basis function information for the
%refinement levels for THB-splines data structure

%INPUT: knotvectors in u,v,w paramteric directions for all the refinement
%levels  - knotvectorU, knotvectorV, knotvectorW
%param- stores the parameters of the registration process in struct format
%plotGrid flag: whether the B-spline grids are to be plotted or not

%OUTPUT: the data structure of the basis array and the element array for
%the hierarchy structure of the B-splines

%Elem:
%'knot_ind','flag_active','IEN','chdElem','cell_centre','node','parElem','actE'
% 1 -- leftmost and bottommost knot index of each element
% 2 -- whther the element is active (1) or not (0)
% 3 -- the non-zero basis functions over the lement at the same level
% 4 -- global indices of the children elements of the current element
% 5 -- centroidal coordinates of the element
% 6 -- nodal indices of the element
% 7 -- parent element index of the current element
% 8 -- active cell index of the element once it is set active

%Basis:
%'basis_ind','flag_active','chdBasis','coeffBasis','suppCell','flag_trunc','flag_ref','flag_bdy','actB'
% 1 -- the i,j,k index of basis function
% 2 -- whther the basis function is active (1) or not (0)
% 3 -- the global indices of the children basis functions
% 4 -- coefficient matrix computed using refinement equation
% 5 -- the elements that have support under the basis function
% 6 -- truncated basis (1) or not (0)
% 7 -- newly inactive basis flag (1) or not activated at all (0)
% 8 -- dirichlet boundary DOF (1) or not (0)
% 9 -- active basis index of the DOF once it is set active

%adding the parameters
maxlevel = param.maxlevel;
nelemU = param.nelemU;
nelemV = param.nelemV;
nelemW = param.nelemW;
nobU = param.nobU;
nobV = param.nobV;
nobW = param.nobW;
Elem = cell(maxlevel,1);
Basis = cell(maxlevel,1);

for level = 1:maxlevel
    fprintf('%d data structure stored\r\n',level);
    
    %Element Data structure
    %knots at level l
    knot_cu = knotvectorU{level,1};
    knot_cv = knotvectorV{level,1};
    knot_cw = knotvectorW{level,1};
    
    %knots at level l+1
    if(level<=maxlevel-1)
        knot_fu = knotvectorU{level+1,1};
        knot_fv = knotvectorV{level+1,1};
        knot_fw = knotvectorW{level+1,1};
    else
        knot_fu = 0;
        knot_fv = 0;
        knot_fw = 0;
    end
    
    %store parent element of each element at level l
    if(level>1)
        Elem{level,7} = parElem;
    else
        Elem{level,7} = -1;
    end
    
    %store element array
    [knot_ind, chdElem, Parcell, IEN, cell_centre, supp_i] = storeElemArray(level,param,knot_cu, knot_cv, knot_cw);
    
    parElem = Parcell;
    
    %knot indices of each element
    Elem{level,1} = knot_ind;
    
    %activity flag of each element: passive or active
    if(level==1)
        Elem{level,2} = ones(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),1);
    else
        Elem{level,2} = zeros(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),1);
    end
    
    %the non zero splines at the current refinement level over each element
    Elem{level,3} = IEN;
    
    %store nodes associated with each element
    Elem{level,6} = -1;
    
    
    %active cell index
    Elem{level,8} = 0;
    
    %store children cells
    if(level<=maxlevel-1)
        Elem{level,4} = chdElem;
    else
        Elem{level,4} = -1*ones(1,8);
    end
    
    %store centroids of each element
    Elem{level,5} = cell_centre;
    
    % Basis Data Structure
    [BB,Chd_Basis,Coeff_Basis,bdyflag] = storeBasisArray(level,param,knot_cu, knot_cv, knot_cw,knot_fu, knot_fv, knot_fw);
    
    %i,j,k index of each B-spline
    Basis{level,1} = BB;
    
    %active flag of the B-spline
    if(level==1)
        Basis{level,2} = ones(nobU(level)*nobV(level)*nobW(level),1);
    else
        Basis{level,2} = zeros(nobU(level)*nobV(level)*nobW(level),1);
    end
    
    %children B-splines
    Basis{level,3} = Chd_Basis;
    
    %coeff matrix of each B-spline
    Basis{level,4} = Coeff_Basis;
    
    %for each B-spline the cell that have the support
    Basis{level,5} = supp_i;
    
    %truncation flag
    Basis{level,6} = zeros(nobU(level)*nobV(level)*nobW(level),1);
    
    %refinement flag
    Basis{level,7} = zeros(nobU(level)*nobV(level)*nobW(level),1);
    %boundary flag
    Basis{level,8} = bdyflag;
    %active B-spline index
    Basis{level,9} = 0;
end
fieldElem = {'knot_ind','flag_active','IEN','chdElem','cell_centre','node','parElem','actE'};
Elem = cell2struct(Elem,fieldElem,2);
fieldBasis = {'basis_ind','flag_active','chdBasis','coeffBasis','suppCell','flag_trunc','flag_ref','flag_bdy','actB'};
Basis = cell2struct(Basis,fieldBasis,2);
disp('finish!');
end