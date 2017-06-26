function [g_phi, g_phiu, g_phiv, g_phiw] = computeGaussPhiMat(u,v,w,param,ders1,ders2,ders3,gg_coeff)

% In this function the basis functions along with the first derivative of
% the basis function in each parametric direction are computed for all the
% Gauss points in each active element

% INPUT:
% u, v, w: the knot indices corresponding to each active element
% param: the struct array containing the parameters
% ders1, ders2, ders3: the values of Nu, Nv, Nw in each parametric
% direction at a particular refinement level
% gg_coeff: the coefficient matrix for the non zero splines in each active
% element

% OUTPUT:
% g_phi: phi = Nu*Nv*Nw computed for all the gauss points in the active
% element
% g_phiu: phiu = Nu'*Nv*Nw computed for all the gauss points in the active
% element
% g_phiv: phiv = Nu*Nv'*Nw computed for all the gauss points in the active
% element
% g_phiw: phiw = Nu*Nv*Nw' computed for all the gauss points in the active
% element
%#codegen

%%
% Degree of each B-spline 
pU = param.pU;
pV = param.pV;
pW = param.pW;

% Gauss order
xlen = param.orderGauss;

% Number of non zero splines over the active element
supp_size = size(gg_coeff,1);


G_phi = zeros(supp_size,xlen,xlen,xlen,'single');
G_phiu = zeros(supp_size,xlen,xlen,xlen,'single');
G_phiv = zeros(supp_size,xlen,xlen,xlen,'single');
G_phiw = zeros(supp_size,xlen,xlen,xlen,'single');

% loop over the gauss points
for gg1 = 1:xlen
    for gg2 = 1:xlen
        for gg3 = 1:xlen         
            
            %compute Nu, Nv, Nw, Nu', Nv', Nw'
            RRD1 = reshape(ders1(u,gg2,:,:),2,(pU+1));
            RRD2 = reshape(ders2(v,gg3,:,:),2,(pV+1));
            RRD3 = reshape(ders3(w,gg1,:,:),2,(pW+1));
            
            %Phi = Nu*Nv*Nw
            RR12 = RRD1(1,:)'*RRD2(1,:);
            RR12 = RR12(:);
            RR123 = RR12*RRD3(1,:);
            RR123 = RR123(:);
            phii = RR123(end:-1:1);
            
            %Phi = Nu'*Nv*Nw
            RRu12 = RRD1(2,:)'*RRD2(1,:);
            RRu12 = RRu12(:);
            RRu123 = RRu12*RRD3(1,:);
            RRu123 = RRu123(:);
            phiiu = RRu123(end:-1:1);
            
            %Phi = Nu*Nv'*Nw
            RRv12 = RRD1(1,:)'*RRD2(2,:);
            RRv12 = RRv12(:);
            RRv123 = RRv12*RRD3(1,:);
            RRv123 = RRv123(:);
            phiiv = RRv123(end:-1:1);
            
            %Phi = Nu*Nv*Nw'
            RRw12 = RRD1(1,:)'*RRD2(1,:);
            RRw12 = RRw12(:);
            RRw123 = RRw12*RRD3(2,:);
            RRw123 = RRw123(:);
            phiiw = RRw123(end:-1:1);
            
            %multiply with the subdivision coefficients
            phi_pi = gg_coeff*phii;
            phi_piu =gg_coeff*phiiu;
            phi_piv = gg_coeff*phiiv;
            phi_piw = gg_coeff*phiiw;

            G_phi(:,gg3,gg2,gg1) = phi_pi;
            G_phiu(:,gg3,gg2,gg1) = phi_piu;
            G_phiv(:,gg3,gg2,gg1) = phi_piv;
            G_phiw(:,gg3,gg2,gg1) = phi_piw;
        end
    end
end

g_phi = G_phi;
g_phiu = G_phiu;
g_phiv = G_phiv;
g_phiw = G_phiw;
end