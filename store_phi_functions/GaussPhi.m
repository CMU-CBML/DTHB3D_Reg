function [PHI,PHIU,PHIV,PHIW,BIGX,BIGY,BIGZ,H] = GaussPhi(ac,Em,knotvectorU,knotvectorV,knotvectorW,Coeff,param)
%#codegen
% This function computes the basis function (phi) along with the first
% derivative in each parametric direction (phi_u, phi_v, phi_w) at each
% gauss point in the active element

% This function also stores the coordinates of the gauss points for the intensities interpolated at the
% gauss points in each active cell stored in BIGX, BIGY, BIGZ

% H matrix stores the (size of each active element/2) at each active element

% INPUT:
% ac: active element array
% Em: element array struct form
% knotvectorU: knot vector in u direction
% knotvectorV: knot vector in v direction
% knotvectorW: knot vector in w direction
% Coeff: coefficient matrix for each active element
% param: the struct variable storing the parameters

% OUTPUT:
% PHI: cell array containing the basis functions values at each gauss point in each active
% element
% PHIU: cell array containing the derivative in u direction
%basis functions values at each gauss point in each active element
% PHIV: cell array containing the derivative in v direction basis functions values at each gauss point in each active
% element
% PHIW: cell array containing the derivative in w direction basis functions values at each gauss point in each active
% element
% BIGX, BIGY, BIGZ: the coordinates of the physical location of the gauss
% points in each active element
% H: the size of each active element divided by 2 (Jacobian of the transformation)
%%
% order of gaussian quadrature
orderGauss = param.orderGauss;

% maximum refinement level
maxlevel = param.maxlevel;

%number of elements at all refinement levels
nelemU = param.nelemU;
nelemV = param.nelemV;
nelemW = param.nelemW;

%degree of splines
pU = param.pU;
pV = param.pV;
pW = param.pW;

%number of splines in each direction
nobU = param.nobU;
nobV = param.nobV;
nobW = param.nobW;
ac_ct = size(ac,1);

%the gauss points of the corresponding gaussian quadrature order
[si1,~] = ggquad(orderGauss);
[si2,~] = ggquad(orderGauss);
[si3,~] = ggquad(orderGauss);

H = zeros(ac_ct,3,'single');

PHI =  cell(ac_ct,1);
PHIU = cell(ac_ct,1);
PHIV = cell(ac_ct,1);
PHIW = cell(ac_ct,1);

DERU = cell(maxlevel,1);
DERV = cell(maxlevel,1);
DERW = cell(maxlevel,1);

for i =1:maxlevel 
    
    
    dersU = zeros(nobU(i),orderGauss,2,pU+1);
    dersV = zeros(nobV(i,1),orderGauss,2,pV+1);
    dersW = zeros(nobW(i,1),orderGauss,2,pW+1);
    
    knotU = knotvectorU{i,1};
    knotV = knotvectorV{i,1};
    knotW = knotvectorW{i,1};
    
    %Compute Nu at the gauss points in u direction
    for j = (pU+1):size(knotU,2)-pU-1
        a1 = knotU(1,j);
        b1 = knotU(1,j+1);
        h1 = (b1-a1)/2;
        m1 = (a1+b1)/2;
        
        s1 = si1.*h1 + m1;
        
        for k=1:size(s1,1)
            index = FindSpan(nobU(i,1)-1,pU,s1(k,1),knotU) + 1;
            if j>size(dersU,1)
                size(dersU)
                pause
            end
            if k>size(dersU,2)
                size(dersU)
                pause
            end
            
            dersU(j,k,:,:) = Der1BasisFun(index-1,s1(k,1),pU,knotU);
        end
    end
    
    %Compute Nv at the gauss points in v direction
    for j = (pV+1):size(knotV,2)-pV-1
        a2 = knotV(1,j);
        b2 = knotV(1,j+1);
        h2 = (b2-a2)/2;
        m2 = (a2+b2)/2;
        
        s2 = si2.*h2 + m2;
        
        for k=1:size(s2,1)
            index = FindSpan(nobV(i,1)-1,pV,s2(k,1),knotV) + 1;
            dersV(j,k,:,:) = Der1BasisFun(index-1,s2(k,1),pV,knotV);
        end
    end
    
    %Compute Nw at the gauss points in w direction
    for j = (pW+1):size(knotW,2)-pW-1
        a3 = knotW(1,j);
        b3 = knotW(1,j+1);
        h3 = (b3-a3)/2;
        m3 = (a3+b3)/2;
        
        s3 = si3.*h3 + m3;
        
        for k=1:size(s3,1)
            index = FindSpan(nobW(i,1)-1,pW,s3(k,1),knotW) + 1;
            dersW(j,k,:,:) = Der1BasisFun(index-1,s3(k,1),pW,knotW);
        end
    end
    
    DERU{i} = dersU;
    DERV{i} = dersV;
    DERW{i} = dersW;
end

BIGX = zeros(orderGauss,ac_ct,orderGauss,orderGauss,'single');
BIGY = zeros(orderGauss,ac_ct,orderGauss,orderGauss,'single');
BIGZ = zeros(orderGauss,ac_ct,orderGauss,orderGauss,'single');

%loop over the active elements
size(Coeff)

parfor i = 1:ac_ct
    
    cell_index = ac(i,1);
    cell_level = ac(i,2);
    
    %coefficient matrix of non zero splines
    gg_coeff = Coeff(i).mat;
    
    %knot vectors at the refinement level of the element
    knotU = knotvectorU{cell_level,1};
    knotV = knotvectorV{cell_level,1};
    knotW = knotvectorW{cell_level,1};
    
    u = Em(cell_level).knot_ind(cell_index,1,:);
    v = Em(cell_level).knot_ind(cell_index,2,:);
    w = Em(cell_level).knot_ind(cell_index,3,:);
    
    a1 = knotU(1,u(1,1));
    b1 = knotU(1,u(1,2));
    h1 = (b1-a1)/2;
    m1 = (a1+b1)/2;
    
    a2 = knotV(1,v(1,1));
    b2 = knotV(1,v(1,2));
    h2 = (b2-a2)/2;
    m2 = (a2+b2)/2;
    
    a3 = knotW(1,w(1,1));
    b3 = knotW(1,w(1,2));
    h3 = (b3-a3)/2;
    m3 = (a3+b3)/2;
    
    %physical coordinate tranformation of the gauss points
    s1 = si1*h1 + m1;
    s2 = si2*h2 + m2;
    s3 = si3*h3 + m3;
    
    [GX,GY,GZ] = meshgrid(s1,s2,s3);
    
    %store the gaussian points in BIGX, BIGY, BIGZ
    % H stores the size of each active element / 2
    BIGX(:,i,:,:) = GX;
    BIGY(:,i,:,:) = GY;
    BIGZ(:,i,:,:) = GZ;
    H(i,:) = [h1,h2,h3];
    
    % Nu, Nv, Nw for the refinement level of the element
    ders1 = DERU{cell_level,1};
    ders2 = DERV{cell_level,1};
    ders3 = DERW{cell_level,1};
    
    %compute phi, phi_u, phi_v, phiw for all the gauss points in the active
    %element
    [g_phi, g_phiu, g_phiv, g_phiw] = computeGaussPhiMat(u(1,1),v(1,1),w(1,1),param,ders1,ders2,ders3,gg_coeff);
    
    %store the phis in the corresponding cell array
    PHI{i} = g_phi;
    PHIU{i} = g_phiu;
    PHIV{i} = g_phiv;
    PHIW{i} = g_phiw;
    
end

BIGX = reshape(BIGX,ac_ct*orderGauss,orderGauss,orderGauss);
BIGY = reshape(BIGY,ac_ct*orderGauss,orderGauss,orderGauss);
BIGZ = reshape(BIGZ,ac_ct*orderGauss,orderGauss,orderGauss);

end