function mainfunction_synthetic_img()
clc;
clear all;
close all;

output_file1 ='output_file_sphere_sun.txt';
fid_out1 = fopen(output_file1,'w');
% if computing Dice Similarity for medical images from Brainweb, setflagDS
% = 1;
setflagDS = 0;
saveVTK = 0;
plotImage = 1; %plot image for display and save the image in .png format
saveImage = 1; %save image result in .mat format

disp('Added paths....');
%add the paths of subfolders of the software
addpaths();
%set the parameters for running registration process
disp('Setting parameters...');
param = setparameters_synthetic();

%open the txt file to enter dice similarity values for each iteration
if(setflagDS ==1)
    output_file='Dice_similarity_sub45_46.txt';
    fid_out = fopen(output_file,'w');
    fprintf(fid_out, 'Level Iteration Background  CSF  GrayMatter  WhiteMatter  Fat  Muscle Muscle/Skin  Skull  vessels  AroundFat  DuraMatter  BoneMarrow\r\n');
end

%tolerance value for stopping criterion in the iteration loop (delta)
tol = 1e-04;
%maximum number of iterations for each refinement level
itermax = 50;
%% 1: Read image data
%Enter file name here
disp('Reading image data...');
load Sphere.mat;
load sun_like.mat;
I1 = Sphere;
I2 = sun_like;
I1_in = I1; %store initial moving image
I2_in = I2; %store initial target image

%Image size
sizeImage = size(I1);

%compute initial RS value
disp('Compute initial RS...');
Idiff = I1-I2;
Idiff2 = (Idiff).^2;
residual_initial = sqrt(sum(sum(sum(Idiff2,3),2),1));
fprintf(fid_out1,'initial residual = %f\r\n',residual_initial);
RS_initial = 1;
fprintf('initial residual %f\n',residual_initial);
fprintf('initial RS value %f\n',RS_initial);

%Store the pixel coordinates
disp('Storing pixel coordinates...');
[pixX, pixY, pixZ]= ndgrid(linspace(1,size(I1,1),size(I1,1)),linspace(1,size(I1,2),size(I1,2)),linspace(1,size(I1,3),size(I1,3)));
X = pixX(:);
Y = pixY(:);
Z = pixZ(:);
pix = [X,Y,Z];

%Initialize the parametric space
disp('Store the knotvectors...');
[CP, knotvectorU, knotvectorV, knotvectorW] = storeKnotArray(param,sizeImage);

%% Compute Elem, Basis, Control point data structure for THBS splines
disp('Set B-spline grid...');

tic
[Dm,Em] = setBsplineGrid_func_withoutNode(knotvectorU,knotvectorV,knotvectorW,param);
toc



%% Store level 1 control points
disp('Compute control points at level 1 using Greville coordinates...');
knotu = knotvectorU{1,1};
knotv = knotvectorV{1,1};
knotw = knotvectorW{1,1};

pp = zeros(param.nobU(1,1)*param.nobV(1,1)*param.nobW(1,1),3);
for k = 1:param.nobW(1,1)
    for j =1:param.nobV(1,1)
        for i = 1:param.nobU(1,1)
            index =(k-1)*param.nobV(1,1)*param.nobU(1,1) + (j-1)*param.nobU(1,1) + i;
            %Greville Abscissae
            coordx = sum(knotu(i+1:i+param.pU))./param.pU;
            coordy = sum(knotv(j+1:j+param.pV))./param.pV;
            coordz = sum(knotw(k+1:k+param.pW))./param.pW;
            pp(index,:) = [coordx,coordy,coordz];
        end
    end
end

%Control points stored for each refienement level
CP(1).pts = pp;
Pm = CP;
Pm_old = Pm;

iterct = 0;
%Loop over each refinement level
disp('Loop over each refinement level...');
for multilev = 0:1:param.maxlevel-1
    [DI1X,DI1Y,DI1Z] = gradient(I1);
    Pm = Pm_old; %Reinitialization of the control points done here
    
    tic
    fprintf('Refinement at level %i...\n',multilev+1);
    if(multilev>0)
        [Em,Dm,Pm] = THB_Refinement(Em,Dm,Pm,knotvectorU, knotvectorV,knotvectorW,bf,CellGrad,meanGrad,param,multilev);
    end
    toc
    
    disp('Collecting active elements, control points and basis functions...');
    [ac, bf, ACP, RHS,Em,Dm] = storeActiveElem(Em,Dm,Pm,multilev);
    ac_ct = size(ac,1);
    bf_ct = size(bf,1);
    fprintf(fid_out1,'active elements = %d, active DOF =  %d at level = %d \r\n',ac_ct,bf_ct,multilev+1);
    %store undeformed control points
    Pm_old = Pm;
    
    disp('Computing the non-zeros spline over each active element and storing coefficient matrices...');
    tic
    [Jm, Coeff] = computeNonZeroSplines_mex(ac, param, Em, Dm);
    toc
    
    disp('Computing the basis functions at pixel coordinates...');
    tic
    numPixels = int64(prod(sizeImage));
    [Pixel, Pix2] = storePixelPhi_mex(numPixels, multilev,pix, knotvectorU, knotvectorV, knotvectorW, Em, Coeff, param);
    for p_ind = 1:numPixels
        %Pixel(p_ind).active_cell=Pix1(p_ind);
        %Pix{p_ind,2} = single(Pix{p_ind,2});
        Pixel(p_ind).phi = Pix2{p_ind};
    end
    clear Pix2
    toc;
    
    tic;
    disp('Computing the basis functions at gaussian points...');
    %compute the gaussian points and weights of the given gauss order
    [Gu,Wu] = ggquad(param.orderGauss);
    [Gv,Wv] = ggquad(param.orderGauss);
    [Gw,Ww] = ggquad(param.orderGauss);
    
    [PHI,PHIU,PHIV,PHIW,BIGX,BIGY,BIGZ,H] = GaussPhi(ac,Em,knotvectorU,knotvectorV,knotvectorW,Coeff,param);
    % interpolate the intensity values of the target image at the gauss
    % points stored in BIGX, BIGY, BIGZ
    cII2 = interp3(pixY, pixX, pixZ, I2, BIGY, BIGX, BIGZ,'*linear',max(I2(:)));
    clear('BIGX','BIGY','BIGZ');
    toc;
    
    %% Start the image registration for each iteration
    disp('Starting the iteration loop for dynamic update of control points...');
    
    %store the initial RHS value in RHS_init
    RHS_init = RHS;
    
    % start the iteration loop
    %% Update the iteration loop here
    
    % Gauss order
    orderGauss = param.orderGauss;
    
    % lambda 1,2 : regularization parameters
    lambda_1 = param.lambda_1;
    lambda_2 = param.lambda_2;
    
    % gamma term in g(x)
    smallNumber = param.smallNumber;
    
    % timestep for each refinement level
    timestep = param.timestep(multilev+1);
    
    % convert the cell array to struct array
    PHI1 = cell2struct(PHI,'mat',2);
    PHIU1 = cell2struct(PHIU,'mat',2);
    PHIV1 = cell2struct(PHIV,'mat',2);
    PHIW1 = cell2struct(PHIW,'mat',2);
    
    iterct_level = 0;
    RS_final = 2;
    
    itermax = 30;
    % while the stopping criterion is satisfied
    while (abs(RS_final-RS_initial)>tol && iterct_level < itermax)
        
        tic
        
        % counting total iterations
        iterct = iterct +1;
        
        % iterations for each refinement level
        iterct_level = iterct_level+1;
        
        RS_initial = RS_final;
        RHS = RHS_init;
        
        % compute the spatial transformation function f(x), f_u(x), f_v(x), f_w(x)
        [BIGXX,BIGYY,BIGZZ,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ] = computenewPoints_mex(Jm,ACP,PHI1,PHIU1,PHIV1,PHIW1,orderGauss);
        
        % interpolate the intensity and grdient at f(x) at the deformed positions of the gauss points
        % cII1:   I1(f(x))
        % cDII1X: I1_x(f(x))
        % cDII1Y: I1_y(f(x))
        % cDII1Z: I1_z(f(x))
        cII1 = interp3(pixY, pixX, pixZ, I1, BIGYY, BIGXX, BIGZZ,'*linear',max(I2(:)));
        cDII1X = interp3(pixY, pixX, pixZ, DI1X,BIGYY, BIGXX, BIGZZ,'*linear',max(I2(:)));
        cDII1Y = interp3(pixY, pixX, pixZ, DI1Y, BIGYY, BIGXX, BIGZZ,'*linear',max(I2(:)));
        cDII1Z = interp3(pixY, pixX, pixZ, DI1Z, BIGYY, BIGXX, BIGZZ,'*linear',max(I2(:)));
        clear('BIGXX','BIGYY','BIGZZ');
        
        % denominator of the fidelity term (g(x))
        denominate = sqrt((cDII1X.^2) + (cDII1Y.^2) + (cDII1Z.^2)+ smallNumber); %g(x)
        
        % fidelity term in x, y, z directions
        Bterm1 = (cII1 - cII2).*2.*cDII1Y./denominate;
        Bterm2 = (cII1 - cII2).*2.*cDII1X./denominate;
        Bterm3 = (cII1 - cII2).*2.*cDII1Z./denominate;
        
        % Now compute the integrals for each basis function in the support of the active cell.
        RHS= compute_Integ_Domain_mex(Jm,Bterm1,Bterm2,Bterm3,BIGMUX,BIGMUY,BIGMUZ,BIGMVX,BIGMVY,BIGMVZ,BIGMWX,BIGMWY,BIGMWZ,RHS,PHI1,PHIU1,PHIV1,PHIW1,lambda_1,lambda_2, Wu,Wv,Ww,H);
        
        % apply dirichlet boundary condition on the control points at the
        % boundary
        RHS = bcondition3D(RHS);
        
        % update the control points P_new = P_old - epsilon*delta(E)
        ACP(:,1:3) = ACP(:,1:3) - timestep.*RHS(:,1:3);
        
        %compute the new positions of the pixel coordiantes using the spatial
        %transformation function
        [pxx, pyy, pzz] = tripleIterLoop_mex(sizeImage, Pixel, Jm, ACP);
        
        %Using the new position of the pixels, compute the new moving image
        I1(I1<0) = 0;
        I1(I1>255) = 255;
        I1 = interp3(pixY, pixX, pixZ, I1,pyy,pxx,pzz,'*linear');
        I1(isnan(I1)) = max(I2(:));
        I1 = round(I1);
        Iplot = interp3(pixY, pixX, pixZ, I1,pyy,pxx,pzz);
        
        %compute the image gradient
        [DI1X,DI1Y,DI1Z] = gradient(I1);
        
        %compute the residual and RS value
        Idiff2 = I1-I2;
        Idiff22 = (Idiff2).^2;
        Sum1 = sum(sum(sum(Idiff22,3),2),1);
        final_residual = sqrt(Sum1);
        RS_final = final_residual/residual_initial;
        fprintf('RS = %f at Iteration = %d   ',RS_final, iterct_level);
        fprintf(fid_out1,'Iteration = %d, residual = %f, RS = %f\r\n',iterct_level, final_residual, RS_final);
        rs(iterct,1) = RS_final;
        iterations(iterct,1) = iterct;
        
        %compute Dice similarity
        if(setflagDS ==1)
            computeDiceSimilarity(fid_out,iterct,I1,I2);
        end
        
        % Plot the image and store the image in .png format
        if(plotImage==1 && saveImage == 1)
            PlotImage(iterct,I1_in,I2_in,Iplot);
        end
        
        if(saveVTK == 1)
            if(mod(iterct,5)==0)
                filename1 = sprintf('../PostProcessing/evolve_image%d.vtk',iterct);
                vtkwrite(filename1,'structured_grid',pixY,pixX,pixZ,'scalars','Intensity',Iplot);
            end
            
        end
        
        if(RS_final>RS_initial)
            break;
        end
        
        toc;
        
    end
    
    %the centroids of the B-spline grid
    cell_co = zeros(ac_ct,3);
    for i = 1:ac_ct
        cell_id = ac(i,1);
        cell_le = ac(i,2);
        cell_co(i,:) = Em(cell_le).cell_centre(cell_id,:);
    end
    
    %compute the absolute difference between the pair of images at centroid
    %of elements
    [CellGrad, meanGrad] = computeDiffGradImage(cell_co,I1,I2, pixX, pixY, pixZ);
    fprintf('Mean Gradient = %f\n',meanGrad);
end

figure;
plot(iterations,rs)
xlabel('Iterations');
ylabel('RS');
end

