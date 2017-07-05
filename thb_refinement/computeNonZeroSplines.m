function[NZ_AElem,Coeff_AElem] = computeNonZeroSplines(ac, param, Em, Dm)

% In this function the non zero splines over each active element from
% different refinement levels are stored. the coefficient matrix of the non
% zero spline over the active elements in stored.

%INPUT:
%ac: active element array
%param: struct array storing the parameters
%Em: element struct array
%Dm: basis function struct array

%OUTPUT:
%NZ_AElem: the non-zero splines stored over each active element
%Coeff_AElem: coefficient matrix of the non zero splines over each active
%element

%%
% degree of B-spline in each direction
pU = param.pU;
pV = param.pV;
pW = param.pW;

% maximum refinement level
maxlevel = param.maxlevel;

%% Non-zero splines for each active cell
% size of the active element array
ac_ct = size(ac,1);
Coeff = cell(ac_ct,1);
Jm = cell(ac_ct,1);

for i = 1:ac_ct %loop over the active elements
    
    counter = 0;
        
    cell_ind = ac(i,1);
    cell_lev = ac(i,2);
    curr_cell = cell_ind;
    
    ctemp = zeros(maxlevel,(pU+1)*(pV+1)*(pW+1),(pU+1)*(pV+1)*(pW+1));
    
    %Collect the non-zeros HB splines over each active cell at a particular
    %refinement level
    local_b = Em(cell_lev,1).IEN(cell_ind,:);
    lb_size = size(local_b,2);
    ct = zeros((pU+1)*(pV+1)*(pW+1));
        
    %count the dimension of cell_sup
    for j = 1:lb_size %loop over non-zero splines at the same level as the element
        if(Dm(cell_lev).flag_active(local_b(1,j),1)==1) %if the spline is active
            counter= counter+1;           
        end
    end
    
     if(cell_lev~=1)
        for m = (cell_lev-1):-1:1 %loop over level in decreasing order till level 1
                        
            curr_cell = Em(m+1).parElem(curr_cell,1); %find the parent cell of the active element
            local_supp = Em(m).IEN(curr_cell,:); %non zero splines of the parent element
            
            ActiveFlagC = Dm(m).flag_active(local_supp,1);                                            
            for j = 1:(pU+1)*(pV+1)*(pW+1)
                %if the spline is active
                if(ActiveFlagC(j)==1)
                    counter=counter+1;                   
                end
            end
        end
     end    
    
    %initialize cell_sup and coef_arr with the correct dimensions
    cell_supp = zeros(counter,1, 'int64');
    coef_arr = zeros(counter,(pU+1)*(pV+1)*(pW+1),'single');
    curr_cell = cell_ind;
    
    counter = 0;
    for j = 1:lb_size
        ct(j,j) = 1;
    end
    ctemp(cell_lev,:,:) = ct;
    
    for j = 1:lb_size %loop over non-zero splines at the same level as the element
        if(Dm(cell_lev).flag_active(local_b(1,j),1)==1) %if the spline is active
            counter= counter+1;
            cell_supp(counter,1) = Dm(cell_lev).actB(local_b(1,j),1); %add to support splines            
            coef_arr(counter,:) = ct(j,:);
        end
    end
    
    
    if(cell_lev~=1)
        for m = (cell_lev-1):-1:1 %loop over level in decreasing order till level 1
            
            local_s1 = Em(m+1).IEN(curr_cell,:);
            curr_cell = Em(m+1).parElem(curr_cell,1); %find the parent cell of the active element
            local_supp = Em(m).IEN(curr_cell,:); %non zero splines of the parent element
            
            cmat = Dm(m).coeffBasis(local_supp,:,:);
            ActiveFlagC = Dm(m).flag_active(local_supp,1);
            ActiveFlagF = Dm(m+1).flag_active(local_s1,1);
            NewInactiveF = Dm(m+1).flag_trunc(local_s1,1);
            
            %compute coefficient matrix of the splines at this level
            [ct,ct11] = computeCoeffMat_mex(pU, pV, pW,ActiveFlagF,NewInactiveF,cmat, local_s1);
            
            ct1 = squeeze(ctemp(m+1,:,:));
            ct = ct*ct1;
            ctnew = ct11*ct1;
            ctemp(m,:,:) = ct;
            
            for j = 1:(pU+1)*(pV+1)*(pW+1)
                %if the spline is active
                if(ActiveFlagC(j)==1)
                    counter=counter+1;
                    %add to support splines
                    cell_supp(counter,1) = Dm(m).actB(local_supp(1,j)); 
          
                    %add the support splines coefficient matrices
                    coef_arr(counter,:) = ctnew(j,:);
                end
            end
        end
    end
    Jm{i,1}= cell_supp; %add the non zero splines at all levels for each active element
    Coeff{i,1} = coef_arr; %add the coefficient matrices of the non zero splines for each active element
end

%convert the cell arrays to struct by hand
NZ_AElem_temp = struct('nzsplines',Jm{1});
NZ_AElem = repmat(NZ_AElem_temp,length(Jm),1);

Coeff_AElem_temp = struct('mat', Coeff{1});
Coeff_AElem = repmat(Coeff_AElem_temp,length(Jm),1);

%Coeff_AElem = struct;
for i=1:length(Jm)
    NZ_AElem(i).nzsplines = Jm{i};
    Coeff_AElem(i).mat = Coeff{i};
end



%NZ_AElem = cell2struct(Jm,'nzsplines',2);
%Coeff_AElem = cell2struct(Coeff,'mat',2);
end