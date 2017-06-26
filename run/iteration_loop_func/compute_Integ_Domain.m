function RHS_final = compute_Integ_Domain(Jm,Bterm1,Bterm2,Bterm3, BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,lambda1,lambda2,w1,w2,w3,H)
%#codegen
% This function computes the integral of the energy functional to update
% the position of the control points

% INPUT:
% Jm: the non-zero splines over the active elements
% Bterm1, Bterm2, Bterm3: the fidelity term in x, y, z direction
% respectively
% BIGMUX, BIGMUY, BIGMUZ: f_u(x) computed at the gauss points in the control
% grid
% BIGMVX, BIGMVY, BIGMVZ: f_v(x) computed at the gauss points in the control
% grid
% BIGMWX, BIGMWY, BIGMWZ: f_w(x) computed at the gauss points in the control
% grid
% RHS: the right hand side of the computation of the update of the control
% points
% PHI1,PHIU1,PHIV1,PHIW1:  the phi and first derivative of the phi in x, y,
% z at the Gauss points in the active elements
% lambda1, lambda2: the regularization parameters
% w1,w2,w3: the weights of the gaussian quadrature
% H: array containing the size of the active elements divided by 2

% OUTPUT:
% RHS_final: update of the right hand side of the equation after the
% computation of delta(E).

%%
ac_ct = size(Jm,1);
bf_ct = size(RHS,1);
xlen = size(Bterm1,2);

parfor i=1:ac_ct
    RHS1 = zeros(bf_ct,4);
    SB = Jm(i).nzsplines;
    supp_phi = PHI1(i).mat;
    supp_phiu = PHIU1(i).mat;
    supp_phiv = PHIV1(i).mat;
    supp_phiw = PHIW1(i).mat;
    hu = H(i,1);
    hv = H(i,2);
    hw = H(i,3);
    
    term1 = Bterm1(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term2 = Bterm2(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term3 = Bterm3(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term4 = BIGMUX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term5 = BIGMUY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term6 = BIGMUZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term7 = BIGMVX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term8 = BIGMVY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term9 = BIGMVZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    term10 = BIGMWX(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term11 = BIGMWY(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    term12 = BIGMWZ(1+(i-1)*xlen:i*xlen,1:xlen,1:xlen);
    
    supp_size = size(supp_phi,1);
    fu_norm = term4.^2 + term5.^2 + term6.^2;
    fv_norm = term7.^2 + term8.^2 + term9.^2;
    fw_norm = term10.^2 + term11.^2 + term12.^2;
    
    fu_fv = term4.*term7+term5.*term8+term6.*term9;
    fv_fw = term7.*term10+term8.*term11+term9.*term12;
    fw_fu = term10.*term4+term11.*term5+term12.*term6;
    
    arx = 2.*(term4.*(fv_norm) - term7.*fu_fv + term4.*fw_norm - term10.*fw_fu);
    brx = 2.*(term7.*(fu_norm) - term4.*fu_fv + term7.*fw_norm - term10.*fv_fw);
    crx = 2.*(term10.*(fv_norm) - term7.*fv_fw + term10.*fu_norm - term4.*fw_fu);
    
    ary = 2.*(term5.*(fv_norm) - term8.*fu_fv + term5.*fw_norm - term11.*fw_fu);
    bry = 2.*(term8.*(fu_norm) - term5.*fu_fv + term8.*fw_norm - term11.*fv_fw);
    cry = 2.*(term11.*(fv_norm) - term8.*fv_fw + term11.*fu_norm - term5.*fw_fu);
    
    arz = 2.*(term6.*(fv_norm) - term9.*fu_fv + term6.*fw_norm - term12.*fw_fu);
    brz = 2.*(term9.*(fu_norm) - term6.*fu_fv + term9.*fw_norm - term12.*fv_fw);
    crz = 2.*(term12.*(fv_norm) - term9.*fv_fw + term12.*fu_norm - term6.*fw_fu);
    
    val1 = zeros(supp_size,1);
    val2 = zeros(supp_size,1);
    val3 = zeros(supp_size,1);
    
    valm1 = zeros(supp_size,xlen,xlen,xlen);
    valm2 = zeros(supp_size,xlen,xlen,xlen);
    valm3 = zeros(supp_size,xlen,xlen,xlen);
    
    for gg1 = 1:xlen
        for gg2 = 1:xlen
            for gg3 = 1:xlen
                
                phi_i  = supp_phi(:,gg1,gg2,gg3);
                phi_ui = supp_phiu(:,gg1,gg2,gg3);
                phi_vi = supp_phiv(:,gg1,gg2,gg3);
                phi_wi = supp_phiw(:,gg1,gg2,gg3);
                
                valm1(:,gg1,gg2,gg3) = (phi_i).*(term1(gg1,gg2,gg3))+ 2*lambda1*(term4(gg1,gg2,gg3).*phi_ui+term7(gg1,gg2,gg3).*phi_vi+term10(gg1,gg2,gg3).*phi_wi) + lambda2*(arx(gg1,gg2,gg3).*phi_ui+brx(gg1,gg2,gg3).*phi_vi+crx(gg1,gg2,gg3).*phi_wi);
                valm2(:,gg1,gg2,gg3) = (phi_i).*(term2(gg1,gg2,gg3))+ 2*lambda1*(term5(gg1,gg2,gg3).*phi_ui+term8(gg1,gg2,gg3).*phi_vi+term11(gg1,gg2,gg3).*phi_wi) + lambda2*(ary(gg1,gg2,gg3).*phi_ui+bry(gg1,gg2,gg3).*phi_vi+cry(gg1,gg2,gg3).*phi_wi);
                valm3(:,gg1,gg2,gg3) = (phi_i).*(term3(gg1,gg2,gg3))+ 2*lambda1*(term6(gg1,gg2,gg3).*phi_ui+term9(gg1,gg2,gg3).*phi_vi+term12(gg1,gg2,gg3).*phi_wi) + lambda2*(arz(gg1,gg2,gg3).*phi_ui+brz(gg1,gg2,gg3).*phi_vi+crz(gg1,gg2,gg3).*phi_wi);
                
                val1(:,1) = val1(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm1(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val2(:,1) = val2(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm2(:,gg1,gg2,gg3).*hu.*hv.*hw;
                val3(:,1) = val3(:,1) + w1(gg1,1).*w2(gg2,1).*w3(gg3,1).*valm3(:,gg1,gg2,gg3).*hu.*hv.*hw;
            end
        end
    end
    
    RHS1(SB,1) =  val1;
    RHS1(SB,2) =  val2;
    RHS1(SB,3) =  val3;
    RHS = RHS + RHS1;
end

RHS_final = RHS;

end