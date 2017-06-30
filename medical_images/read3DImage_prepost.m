function [Moving,Fixed] = read3DImage_prepost(filenameI0,filenameI1)

%This function reads the .mha files files of brain MR data taken from SPL lab data base (http://www.spl.harvard.edu/publications/item/view/1915)
% to .mat files as input for the registration process
%INPUT: filenameI0: the filename of the source image
%filenameI1: the filename of the target image
%OUTPUT: Moving: matrix of the source image
%Fixed: matrix of the target image
% Image is padded with zeros on the boundaries to make the interpolation
% exact


%read raw binary files
img1 = ReadData3D(filenameI0);
img2 = ReadData3D(filenameI1);
[BrainI0] = double(img1);
[BrainI1] = double(img2);
BrainI0 = resize3Dmatrix(size(BrainI1,1),size(BrainI1,2),size(BrainI1,3),BrainI0);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %convert image to 0-255
BrainI0=(BrainI0-min(BrainI0(:)))/(max(BrainI0(:))-min(BrainI0(:)))*255;
BrainI1=(BrainI1-min(BrainI1(:)))/(max(BrainI1(:))-min(BrainI1(:)))*255;
%
Moving = permute(BrainI0,[1,3,2]); %moving Image
Fixed = permute(BrainI1,[1,3,2]);  %fixed Image
end