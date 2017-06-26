function [BB,Chd_Basis,Coeff_Basis,bdyflag] = storeBasisArray(level,param,knot_cu, knot_cv, knot_cw,knot_fu, knot_fv, knot_fw)

% This function stores the parts of the element struct array for a
% particular refinement level
%INPUT:
% level = current refinement level
% param = parametric struct array
% knot_cu = knotvector of the current refinement level in u direction
% knot_cv = knotvector of the current refinement level in v direction
% knot_cw = knotvector of the current refinement level in w direction
% knot_fu = knotvector of the next refinement level in u direction
% knot_fv = knotvector of the next refinement level in v direction
% knot_fw = knotvector of the next refinement level in w direction
%OUTPUT:
% BB =  indices array for the DOF at the current refinement level
% Chd_basis = children basis function indices for the splines at the current refinement
%level
% Coeff_Basis = coefficient matrix using the refinement equation of the splines at current refinement level
% bdyflag = DOF satisfying the Dirichlet boundary condition

bc = 0;
maxlevel = param.maxlevel;
nobU = param.nobU;
nobV = param.nobV;
nobW = param.nobW;
pU = param.pU;
pV = param.pV;
pW = param.pW;

bdyflag = zeros(nobU(level)*nobV(level)*nobW(level),1);

BB = zeros(nobU(level)*nobV(level)*nobW(level),3);

if(level<maxlevel)
    Chd_Basis = zeros(nobU(level)*nobV(level)*nobW(level),(pU+2)*(pV+2)*(pW+2));
    Coeff_Basis = zeros(nobU(level)*nobV(level)*nobW(level),(pU+2)*(pV+2)*(pW+2),2);
else
    Chd_Basis  = -1*ones(1,(pU+2)*(pV+2)*(pW+2));
    Coeff_Basis = -1*ones(1,(pU+2)*(pV+2)*(pW+2),2);
end

for basis_k =1:nobW(level,1)
    for basis_j = 1:nobV(level,1)
        for basis_i = 1:nobU(level,1)
            bc=bc+1;
            
            BB(bc,:) = [basis_i,basis_j,basis_k];
            
            if(basis_i==1 || basis_i==nobU(level,1) || basis_j==1 || basis_j==nobV(level,1)|| basis_k==1 || basis_k==nobW(level,1)),
                bdyflag(bc,1) = 1;
            end
            
            if(level<=maxlevel-1)
                intc_u1 = knot_cu(1,basis_i);
                intc_u2 = knot_cu(1,basis_i+pU+1);
                
                intc_v1 = knot_cv(1,basis_j);
                intc_v2 = knot_cv(1,basis_j+pV+1);
                
                intc_w1 = knot_cw(1,basis_k);
                intc_w2 = knot_cw(1,basis_k+pW+1);
                
                if(basis_i < pU+1)
                    uindex_start = basis_i;
                else
                    uindex_start = FindSpan(nobU(level+1,1),pU,intc_u1,knot_fu)+1;
                end
                
                if((nobU(level,1)-basis_i)<pU)
                    uindex_end = size(knot_fu,2)-(nobU(level,1)-basis_i);
                else
                    uindex_end = FindSpan(nobU(level+1,1),pU,intc_u2,knot_fu)+1;
                end
                
                if(basis_j < pV+1)
                    vindex_start = basis_j;
                else
                    vindex_start = FindSpan(nobV(level+1,1),pV,intc_v1,knot_fv)+1;
                end
                
                if((nobV(level,1)-basis_j)<pV)
                    vindex_end = size(knot_fv,2)-(nobV(level,1)-basis_j);
                else
                    vindex_end = FindSpan(nobV(level+1,1),pV,intc_v2,knot_fv)+1;
                end
                
                if(basis_k < pW+1)
                    windex_start = basis_k;
                else
                    windex_start = FindSpan(nobW(level+1,1),pW,intc_w1,knot_fw)+1;
                end
                
                if((nobW(level,1)-basis_k)<pW)
                    windex_end = size(knot_fw,2)-(nobW(level,1)-basis_k);
                else
                    windex_end = FindSpan(nobW(level+1,1),pW,intc_w2,knot_fw)+1;
                end
                
                Knotu = knot_cu(basis_i:(basis_i+pU+1));
                Knotv = knot_cv(basis_j:(basis_j+pV+1));
                Knotw = knot_cw(basis_k:(basis_k+pW+1));
                
                newKnotu = knot_fu(uindex_start:uindex_end);
                newKnotv = knot_fv(vindex_start:vindex_end);
                newKnotw = knot_fw(windex_start:windex_end);
                
                nob_childu = size(newKnotu,2)-pU-1;
                nob_childv = size(newKnotv,2)-pV-1;
                nob_childw = size(newKnotw,2)-pW-1;
                
                nu = size(Knotu,2)-pU-1;
                nv = size(Knotv,2)-pV-1;
                nw = size(Knotw,2)-pW-1;
                
                TmatU =  Tmatrix(Knotu,newKnotu,pU);
                TmatV =  Tmatrix(Knotv,newKnotv,pV);
                TmatW =  Tmatrix(Knotw,newKnotw,pW);
                
                TmatU = TmatU(1:nob_childu,1:nu);
                TmatV = TmatV(1:nob_childv,1:nv);
                TmatW = TmatW(1:nob_childw,1:nw);
                
                cc=0;
                
                for c_childw = 1:nob_childw
                    for c_childv = 1:nob_childv
                        for c_childu = 1:nob_childu
                            
                            cc= cc+1;
                            cc1 = uindex_start+(c_childu-1);
                            cc2 = vindex_start+(c_childv-1);
                            cc3 = windex_start+(c_childw-1);
                            
                            cbb = nobU(level+1,1)*nobV(level+1,1)*(cc3-1) + nobU(level+1,1)*(cc2-1) + cc1;
                            Chd_Basis(bc,cc) = cbb;
                            Coeff_Basis(bc,cc,1) = TmatU(c_childu,1)*TmatV(c_childv,1)*TmatW(c_childw,1);
                            Coeff_Basis(bc,cc,2) = cbb;
                            
                        end
                    end
                end
            end
        end
    end
end
end