function [Elem_final,Basis_final,P_final] = THB_Refinement(Elem,Basis,Cpt,knotvectorU, knotvectorV,knotvectorW,ActiveBasis,CellGrad,meanGrad,param,multilev)

%convert struct aray to cell array
% the reason is refinement function is faster when it works on cell arrays
% than struct data structure
Em = struct2cell(Elem)';
Dm = struct2cell(Basis)';
Pm = struct2cell(Cpt)';
bf_ct = size(ActiveBasis,1);
bf = ActiveBasis;

%degree of B-splines
pU = param.pU;
pV = param.pV;
pW = param.pW;

%loop over the active B-splines
for j =1:bf_ct
    
    %B-spline index, B-spline level
    bbc = bf(j,1);
    bf_lev = bf(j,2);
    
    %refinement parameter
    rho = param.rho(multilev);
    
    %loop over the support cells of the active splines
    supp_cells = Basis(bf_lev).suppCell(bbc,:);
    
    %compute the average image difference value over the supporting
    %elements of the active B-spline
    grad = 0;
    supp_ct = 0;
    for i =1:size(supp_cells,2)
        if(supp_cells(1,i)~=0)
            supp_ct = supp_ct + 1;
            ac_ind = Elem(bf_lev).actE(supp_cells(1,i),1);
            grad  = grad + CellGrad(ac_ind,1);
        end
    end
    grad = grad/supp_ct;
    
    %Refinement to create next level
    %if the average value of Ig over the support elements is greater than
    %mean image difference times rho, REFINE spline
    if(grad>=rho*meanGrad)
        [Dm,Em,Pm] =  Refine3D_truncation(bbc,bf_lev,Dm,Em,Pm,knotvectorU,knotvectorV,knotvectorW,pU,pV,pW);
    end
end

%convert cell array to struct array 
fieldElem = {'knot_ind','flag_active','IEN','chdElem','cell_centre','node','parElem','actE'};
Elem_final = cell2struct(Em,fieldElem,2);
fieldBasis = {'basis_ind','flag_active','chdBasis','coeffBasis','suppCell','flag_trunc','flag_ref','flag_bdy','actB'};
Basis_final = cell2struct(Dm,fieldBasis,2);
P_final = cell2struct(Pm,'pts',2);
end