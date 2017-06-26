function [Pixel, Pix2] = storePixelPhi(numPixels, multilev,pixel, knotvectorU, knotvectorV, knotvectorW, Em, Coeff, param)
%#coder
% this function computes the basis function values (phi) at the coordinates
% of the pixels

% INPUT:
% I1: moving image
% multilev: refinement level -1
% pixel: coordinates of the pixels [x,y,z]
% knotvectorU: knot vectors in u direction
% knotvectorV: knot vectors in v direction
% knotvectorW: knot vectors in w direction
% Em: the element array struct variable
% Coeff: coefficient matrix of non zero spline sover each active element
% param: the struct variable of all the parameters

% OUTPUT:
% Pix: cell array of size = number of pixels
% [1]: active cell index
% [2]: phi matrix, values of basis functions of the non zero splines

%% Store phis at pixel coordinates
%Pix1 = zeros(numPixels,1);
Pix2 = cell(numPixels,1);

%knot vectors
kuMAX = knotvectorU{multilev+1,1};
kvMAX = knotvectorV{multilev+1,1};
kwMAX = knotvectorW{multilev+1,1};

%number of basis functions
nobU = param.nobU;
nobV = param.nobV;
nobW = param.nobW;

% degree of splines
pU = param.pU;
pV = param.pV;
pW = param.pW;
nelemU = param.nelemU;
nelemV = param.nelemV;

% number of basis functions of refinement level = l
unobMAX = nobU(multilev+1,1);
vnobMAX = nobV(multilev+1,1);
wnobMAX = nobW(multilev+1,1);


% loop over the pixels
Pixel_temp = struct('active_cell',0,'phi',single([]));
Pixel =repmat(Pixel_temp,numPixels,1);
%Pixel(numPixels) = struct('active_cell',0,'phi',single([]));

Pix2 = coder.nullcopy(Pix2);
parfor p_ind = 1:numPixels
    [act_ind, phi_pi]=storePixelPhiBody(p_ind, multilev,pixel, knotvectorU, knotvectorV, knotvectorW, kuMAX, kvMAX, kwMAX, pU, pV, pW, nelemU, nelemV, nobU, nobV, nobW, unobMAX, vnobMAX, wnobMAX, Em, Coeff);
    %store the active element and the corresponding phi values
    Pixel(p_ind).active_cell = act_ind;
    Pix2{p_ind} = phi_pi;
end
end

function [act_ind, phi_pi]=storePixelPhiBody(p_ind, multilev,pixel, knotvectorU, knotvectorV, knotvectorW, kuMAX, kvMAX, kwMAX, pU, pV, pW, nelemU, nelemV, nobU, nobV, nobW, unobMAX, vnobMAX, wnobMAX, Em, Coeff)
%find the knot index corresponding to the pixel coordinate
uu = FindSpan(unobMAX-1, pU, pixel(p_ind,1), kuMAX) + 1;
vv = FindSpan(vnobMAX-1, pV, pixel(p_ind,2), kvMAX) + 1;
ww = FindSpan(wnobMAX-1, pW, pixel(p_ind,3), kwMAX) + 1;

%find the cell index correspoding to the knot index
cellx = (uu-pU);
celly = (vv-pV);
cellz = (ww-pW);

cell_ind = nelemU(multilev+1,1)*nelemV(multilev+1,1)*(cellz-1)+nelemU(multilev+1,1)*(celly-1)+cellx;
curr_cell = cell_ind;

%if the cell is active at that refinement level
act_ind = 0;
knotuu = kuMAX;
knotvv = kvMAX;
knotww = kwMAX;
if(Em(multilev+1,1).flag_active(cell_ind,1) == 1)
    %store the active cell index
    act_ind = Em(multilev+1,1).actE(cell_ind,1);
    knotuu = kuMAX;
    knotvv = kvMAX;
    knotww = kwMAX;
    
else
    
    for m = (multilev+1):-1:2
        % go to the parent cell of the inactive cell, check if it is
        % active
        curr_cell = Em(m,1).parElem(curr_cell,1);
        
        if(Em((m-1),1).flag_active(curr_cell,1) == 1)
            %store the active cell index
            act_ind = Em((m-1),1).actE(curr_cell,1);
            knotuu = knotvectorU{(m-1),1};
            knotvv = knotvectorV{(m-1),1};
            knotww = knotvectorW{(m-1),1};
            
            unob1 = nobU(m-1,1);
            vnob1 = nobV(m-1,1);
            wnob1 = nobW(m-1,1);
            
            uu = FindSpan(unob1-1,pU,pixel(p_ind,1),knotuu) + 1;
            vv = FindSpan(vnob1-1,pV,pixel(p_ind,2),knotvv) + 1;
            ww = FindSpan(wnob1-1,pW,pixel(p_ind,3),knotww) + 1;
            break;
        end
        
    end
end

%coeeficient matrix of the activ element selected for the pixel
%coordinate
pix_coeff = Coeff(act_ind).mat;

%compute the basis function
RR1 = BasisFun((uu-1),pixel(p_ind,1),pU,knotuu);
RR2 = BasisFun((vv-1),pixel(p_ind,2),pV,knotvv);
RR3 = BasisFun((ww-1),pixel(p_ind,3),pW,knotww);

RR12 = RR1(1,:)'*RR2(1,:);
RR12 = RR12(:);
RR123 = RR12*RR3(1,:);
RR123 = RR123(:);
phii = RR123(end:-1:1);

% compute phi = coeff*Nu*Nv*Nw
phi_pi = pix_coeff*phii;
end
