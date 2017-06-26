function [BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ] = computenewPoints(Jm,ACP,PHI1,PHIU1,PHIV1,PHIW1,xlen)
%#codegen
% This function generates the coordinates of the gauss points in the
% deformed control grid i.e. f(x), f_u(x), f_v(x) and f_w(x)

% INPUT:
% Jm,ACP,PHI1,PHIU1,PHIV1,PHIW1,xlen
% Jm: the non zero splines over each active element
% ACP: the active control points
% PHI1: the basis function values at the gauss points in each active
% element
% PHIU1,PHIV1,PHIW1: the first derivative of the basis function at each
% gauss points in each active element
% xlen: the gauss quadrature order

% OUTPUT:
% BIGXX, BIGYY, BIGZZ: f(x) computed at the gauss points in the control
% grid
% BIGMUX, BIGMUY, BIGMUZ: f_u(x) computed at the gauss points in the control
% grid
% BIGMVX, BIGMVY, BIGMVZ: f_v(x) computed at the gauss points in the control
% grid
% BIGMWX, BIGMWY, BIGMWZ: f_w(x) computed at the gauss points in the control
% grid

%% 
%number of active elements
ac_ct = length(Jm);

BIGXX = zeros(xlen,ac_ct,xlen,xlen,'single');
BIGYY = zeros(xlen,ac_ct,xlen,xlen,'single');
BIGZZ = zeros(xlen,ac_ct,xlen,xlen,'single');

BIGMUX = zeros(xlen,ac_ct,xlen,xlen,'single');
BIGMUY = zeros(xlen,ac_ct,xlen,xlen,'single');
BIGMUZ = zeros(xlen,ac_ct,xlen,xlen,'single');

BIGMVX = zeros(xlen,ac_ct,xlen,xlen,'single');
BIGMVY = zeros(xlen,ac_ct,xlen,xlen,'single');
BIGMVZ = zeros(xlen,ac_ct,xlen,xlen,'single');

BIGMWX = zeros(xlen,ac_ct,xlen,xlen,'single');
BIGMWY = zeros(xlen,ac_ct,xlen,xlen,'single');
BIGMWZ = zeros(xlen,ac_ct,xlen,xlen,'single');

%loop over active elements
parfor i = 1:ac_ct
    % the non zero splines over each active element
    SB = Jm(i).nzsplines;
    supp_size = size(SB,1);
    %control points of the non zero splines
    pts = ACP(SB,1:3);
    
    %phi, phiu, phiv, phiw of the gauss points
    supp_phi = PHI1(i).mat;
    supp_phiu = PHIU1(i).mat;
    supp_phiv = PHIV1(i).mat;
    supp_phiw = PHIW1(i).mat;
    
    supp_phi = reshape(supp_phi,supp_size,xlen*xlen*xlen);
    supp_phiu = reshape(supp_phiu,supp_size,xlen*xlen*xlen);
    supp_phiv = reshape(supp_phiv,supp_size,xlen*xlen*xlen);
    supp_phiw = reshape(supp_phiw,supp_size,xlen*xlen*xlen);
    
    %Now to compute the vectors, fx, fux, fvx
    SBX = pts(:,1)'*supp_phi;
    SBY = pts(:,2)'*supp_phi;
    SBZ = pts(:,3)'*supp_phi;
    
    SBUX = pts(:,1)'*supp_phiu;
    SBUY = pts(:,2)'*supp_phiu;
    SBUZ = pts(:,3)'*supp_phiu;
    
    SBVX = pts(:,1)'*supp_phiv;
    SBVY = pts(:,2)'*supp_phiv;
    SBVZ = pts(:,3)'*supp_phiv;
    
    SBWX = pts(:,1)'*supp_phiw;
    SBWY = pts(:,2)'*supp_phiw;
    SBWZ = pts(:,3)'*supp_phiw;
    
    SBX = reshape(SBX,xlen,xlen,xlen);
    SBY = reshape(SBY,xlen,xlen,xlen);
    SBZ = reshape(SBZ,xlen,xlen,xlen);
    
    SBUX = reshape(SBUX,xlen,xlen,xlen);
    SBUY = reshape(SBUY,xlen,xlen,xlen);
    SBUZ = reshape(SBUZ,xlen,xlen,xlen);
    
    SBVX = reshape(SBVX,xlen,xlen,xlen);
    SBVY = reshape(SBVY,xlen,xlen,xlen);
    SBVZ = reshape(SBVZ,xlen,xlen,xlen);
    
    SBWX = reshape(SBWX,xlen,xlen,xlen);
    SBWY = reshape(SBWY,xlen,xlen,xlen);
    SBWZ = reshape(SBWZ,xlen,xlen,xlen);
    
    BIGXX(:,i,:,:) = SBX;
    BIGYY(:,i,:,:) = SBY;
    BIGZZ(:,i,:,:) = SBZ;
    
    BIGMUX(:,i,:,:) = SBUX;
    BIGMUY(:,i,:,:) = SBUY;
    BIGMUZ(:,i,:,:) = SBUZ;
    
    BIGMVX(:,i,:,:) = SBVX;
    BIGMVY(:,i,:,:) = SBVY;
    BIGMVZ(:,i,:,:) = SBVZ;
    
    BIGMWX(:,i,:,:) = SBWX;
    BIGMWY(:,i,:,:) = SBWY;
    BIGMWZ(:,i,:,:) = SBWZ;
end

BIGXX = reshape(BIGXX,ac_ct*xlen,xlen,xlen);
BIGYY = reshape(BIGYY,ac_ct*xlen,xlen,xlen);
BIGZZ = reshape(BIGZZ,ac_ct*xlen,xlen,xlen);

BIGMUX = reshape(BIGMUX,ac_ct*xlen,xlen,xlen);
BIGMUY = reshape(BIGMUY,ac_ct*xlen,xlen,xlen);
BIGMUZ = reshape(BIGMUZ,ac_ct*xlen,xlen,xlen);

BIGMVX = reshape(BIGMVX,ac_ct*xlen,xlen,xlen);
BIGMVY = reshape(BIGMVY,ac_ct*xlen,xlen,xlen);
BIGMVZ = reshape(BIGMVZ,ac_ct*xlen,xlen,xlen);

BIGMWX = reshape(BIGMWX,ac_ct*xlen,xlen,xlen);
BIGMWY = reshape(BIGMWY,ac_ct*xlen,xlen,xlen);
BIGMWZ = reshape(BIGMWZ,ac_ct*xlen,xlen,xlen);