function [AElem, ABasis,ACpts, ARHS,Elem,Basis] = storeActiveElem(Elem,Basis,P,multilev)
% This function stores the active basis functions, elements and control
% points after refinement. It also stores the active indices for each
% element in the array structure which is set active

%INPUT:
%Elem: element array structure
%Basis: basis function array structure
%P: control points array structure
%multilev: current refinement level -1 

%OUPUT:
%AElem: active element array
%ABasis: active basis function array
%ACpts: active control points
%ARHS: right hand side for the control point update
%Elem1: element array structure after storing the active cell indices
%Basis1: basis array structure after storing the active basis function
%indices

%#codegen
ac_ct = 0;
bf_ct = 0;
ACP = [];
ac = [];
bf = [];

for level = 1:(multilev+1)
    %load the active flags for the elements and B-splines at a particular level 
    sizee = size(Elem(level).flag_active,1);
    sizeb = size(Basis(level).flag_active,1);
    %load the control points at a particular level
    PP = P(level).pts;
    
    for i = 1:sizee
        %if the element is set as active (1)
        if(Elem(level).flag_active(i,1)==1)
            ac_ct = ac_ct + 1;
            %add the element index and level to the active cell array
            ac = [ac;i,level];
            %update the active cell index for the element in hierarchy
            %array structure
            Elem(level).actE(i,1) = ac_ct;
        end
    end
    
    for j = 1:sizeb
        %if the basis function is set as active (1)
        if(Basis(level).flag_active(j,1)==1)
            bf_ct = bf_ct+1;
            %add the B-spline index and level to the active basis array
            bf = [bf;j,level];
            %update the active basis function index for the basis in hierarchy
            %array structure
            Basis(level).actB(j,1) = bf_ct;
            %add the control point coordinates and boundary flag to the
            %active control point array
            ACP = [ACP;PP(j,:),Basis(level).flag_bdy(j,1)];
        end
    end
end

bf_ct = size(bf,1);
AElem = ac;
ABasis = bf;
ACpts = ACP;
ARHS = zeros(bf_ct,4);
ARHS(:,4) = ACP(:,4);
end