function [cK, chdE,PC, BBv, Cellc,supcell_basis] = storeElemArray(level,param,knot_cu, knot_cv, knot_cw)
% This function stores the parts of the element struct array for a
% particular refinement level
%INPUT:
% level = current refinement level
% param = parametric struct array
% knot_cu = knotvector of the current refinement level in u direction
% knot_cv = knotvector of the current refinement level in v direction
% knot_cw = knotvector of the current refinement level in w direction
%OUTPUT:
% cK = knot indices array for the elements at the current refinement level
% chdE = children element indices for the elements at the current refinement
%level
% PC = parent element index of the elements at current refinement level
% BBv = the non-zero DOF for the elements at current refinement level
% Cellc = centroids of the elements at the current refinement level
% suppcell_basis = the elements under the support of each basis function

%parameters
maxlevel = param.maxlevel;

nobU = param.nobU;
nobV = param.nobV;
nobW = param.nobW;

nelemU = param.nelemU;
nelemV = param.nelemV;
nelemW = param.nelemW;

pU = param.pU;
pV = param.pV;
pW = param.pW;

kV = param.kV;
kW = param.kW;


supp_ct = zeros(nobU(level,1)*nobV(level,1)*nobW(level,1),1);
supp_i  = zeros(nobU(level,1)*nobV(level,1)*nobW(level,1),(1+pU)*(1+pV)*(1+pW));

glnumnodes = 0;
E = 0;

if(level<=maxlevel-1)
    Parcell = zeros(nelemU(level+1,1)*nelemV(level+1,1)*nelemW(level+1,1),1);
    PC = zeros(nelemU(level+1,1)*nelemV(level+1,1)*nelemW(level+1,1),1);
else
    Parcell = zeros(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),1);
    PC = zeros(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),1);
end

BBvector = zeros(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),(1+pU)*(1+pV)*(1+pW));
connectKnots = zeros(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),3,2);
Cell_centre = zeros(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),3);
chdElem= zeros(nelemU(level,1)*nelemV(level,1)*nelemW(level,1),8);


for intW = 1:nobW(level,1)
    for intV = 1:nobV(level,1)
        for intU = 1:nobU(level,1)
            
            glnumnodes = glnumnodes + 1;
            
            if((intU>=(pU+1))&&(intV>=(pV+1))&&(intW>=(pW+1)))
                
                E = E+1;
                i = 1;
                
                connectKnots(E,1,:) = [intU, intU+1];
                connectKnots(E,2,:) = [intV, intV+1];
                connectKnots(E,3,:) = [intW, intW+1];
                
                xind = intU-pU;
                yind = intV-pV;
                zind = intW-pW;
                
                if(level<=maxlevel-1)
                    
                    chdElem(E,1) = nelemU(level+1,1)*nelemV(level+1,1)*((2*zind-1)-1)+nelemU(level+1,1)*((2*yind-1)-1)+ (2*xind-1);
                    chdElem(E,2) = nelemU(level+1,1)*nelemV(level+1,1)*((2*zind-1)-1)+nelemU(level+1,1)*((2*yind-1)-1)+ (2*xind);
                    chdElem(E,3) = nelemU(level+1,1)*nelemV(level+1,1)*((2*zind-1)-1)+nelemU(level+1,1)*(2*yind-1)+ (2*xind-1);
                    chdElem(E,4) = nelemU(level+1,1)*nelemV(level+1,1)*((2*zind-1)-1)+nelemU(level+1,1)*(2*yind-1)+ (2*xind);
                    
                    chdElem(E,5) = nelemU(level+1,1)*nelemV(level+1,1)*(2*zind-1)+nelemU(level+1,1)*((2*yind-1)-1)+ (2*xind-1);
                    chdElem(E,6) = nelemU(level+1,1)*nelemV(level+1,1)*(2*zind-1)+nelemU(level+1,1)*((2*yind-1)-1)+ (2*xind);
                    chdElem(E,7) = nelemU(level+1,1)*nelemV(level+1,1)*(2*zind-1)+nelemU(level+1,1)*(2*yind-1)+ (2*xind-1);
                    chdElem(E,8) = nelemU(level+1,1)*nelemV(level+1,1)*(2*zind-1)+nelemU(level+1,1)*(2*yind-1)+ (2*xind);
                    
                    for pct = 1:8
                        Parcell(chdElem(E,pct),1) = E;
                    end
                    
                end
                
                for loci = 1:pU+1
                    for locj = 1:pV+1
                        for lock = 1:pW+1
                            B  = glnumnodes - (loci-1)*(kW(level,1)-pW-1)*(kV(level,1)-pV-1)-(locj-1)*(kW(level,1)-pW-1) -lock + 1;
                            BBvector(E,i) = B;
                            
                            supp_ct(B,1) = supp_ct(B,1)+1;
                            supp_i(B,supp_ct(B,1)) = E;
                            i=i+1;
                        end
                    end
                end
                
                Cell_centre(E,1) = knot_cu(intU)+(knot_cu(intU+1)-knot_cu(intU))/2;
                Cell_centre(E,2)= knot_cv(intV)+(knot_cv(intV+1)-knot_cv(intV))/2;
                Cell_centre(E,3) = knot_cw(intW)+(knot_cw(intW+1)-knot_cw(intW))/2;
                
            end
        end
    end
end

cK = connectKnots;
chdE = chdElem;
PC = Parcell;
BBv = BBvector;
Cellc = Cell_centre;
supcell_basis = supp_i;
end