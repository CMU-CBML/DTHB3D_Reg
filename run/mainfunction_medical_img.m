function mainfunction_medical_img()
clear all;
close all;

% if computing Dice Similarity for medical images from Brainweb, setflagDS
% = 1;
setflagDS = 1;
plotImage = 1; %plot image for display and save the image in .png format
saveImage = 1; %save image result in .mat format

disp('Added paths....');
%add the paths of subfolders of the software
addpaths();
%set the parameters for running registration process
disp('Setting parameters...');
param = setparameters_brainweb();

%open the txt file to enter dice similarity values for each iteration
if(setflagDS ==1)
    output_file='Dice_similarity_sub45_46.txt';
    fid_out = fopen(output_file,'w');
    fprintf(fid_out, 'Level Iteration Background  CSF  GrayMatter  WhiteMatter  Fat  Muscle Muscle/Skin  Skull  vessels  AroundFat  DuraMatter  BoneMarrow\r\n');
end

output_file1 ='output_file_brainweb.txt';
fid_out1 = fopen(output_file1,'w');
%tolerance value for stopping criterion in the iteration loop (delta)
tol = 1e-06;
%maximum number of iterations for each refinement level
itermax = 50;

%% 1: Read image data
%Enter file name here
disp('Reading image data...');
filenameI0='subject45_crisp_v.rawb';
filenameI1='subject46_crisp_v.rawb';

[I1,I2] = read3DImage_brainWeb(filenameI0,filenameI1);
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
disp('Compute control points at level 1 using Greville abscissae...');
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
rs = zeros(1000,1);
iterations = zeros(1000,1);

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
    tic
    [ac, bf, ACP, RHS,Em,Dm] = storeActiveElem(Em,Dm,Pm,multilev);
    ac_ct = size(ac,1);
    bf_ct = size(bf,1);
    fprintf(fid_out1,'active elements = %d, active DOF =  %d at level = %d \r\n',ac_ct,bf_ct,multilev+1);
    toc
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
    cII2 = interp3(pixY, pixX, pixZ, I2, BIGY, BIGX, BIGZ,'*spline');
    clear('BIGX','BIGY','BIGZ');
    toc;
    
    %% Start the image registration for each iteration
    disp('Starting the iteration loop for dynamic update of control points...');
    
    %store the initial RHS value in RHS_init
    RHS_init = RHS;
    
    % start the iteration loop
    IterationLoop
    
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
plot(rs,iterations)
end

