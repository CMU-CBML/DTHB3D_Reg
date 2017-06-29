function [Moving,Fixed] = read3DImage_brainWeb(filenameI0,filenameI1)
%source code taken from the website https://github.com/stellaccl/cdmffd-image-registration
%Paper: Chan, Chiu Ling, et al. "Two and Three Dimensional Image Registration Based on B-Spline Composition
% and Level Sets." Communications in Computational Physics 21.2 (2017): 600-622.
%This function reads the raw binary brain MRI files taken from Brainweb
%wesite and converts to .mat files as input for the registration process 
%INPUT: filenameI0: the filename of the source image 
%filenameI1: the filename of the target image
%OUTPUT: Moving: matrix of the source image
%Fixed: matrix of the target image
% Image is padded with zeros on the boundaries to make the interpolation
% exact

Moving = zeros(200,224,200);
Fixed = zeros(200,224,200);

%read raw binary files
[BrainI0] = readrawb(filenameI0);
[BrainI1] = readrawb(filenameI1);

scale=1/2;
nx=size(BrainI0,2);
ny=size(BrainI0,1);
nz=size(BrainI0,3);
new_size_x=round(nx*scale);
new_size_y=round(ny*scale);
new_size_z=round(nz*scale);

BrainI0 = resize3Dmatrix(new_size_x,new_size_y,new_size_z,BrainI0);
BrainI1 = resize3Dmatrix(new_size_x,new_size_y,new_size_z,BrainI1);

BrainI0=(BrainI0-min(BrainI0(:)))/(max(BrainI0(:))-min(BrainI0(:)))*11;
BrainI0=round(BrainI0);

BrainI1=(BrainI1-min(BrainI1(:)))/(max(BrainI1(:))-min(BrainI1(:)))*11;
BrainI1=round(BrainI1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert image to 0-255
Img1=BrainI0*23.1818;
Img2=BrainI1*23.1818;

Moving(11:191,4:220,11:191) = Img1; %moving Image
Fixed(11:191,4:220,11:191) = Img2;  %fixed Image
end