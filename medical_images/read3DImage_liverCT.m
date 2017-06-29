function [Moving,Fixed] = read3DImage_liverCT(filenameI0,filenameI1)
%This function reads the mha image volume liver CT files taken from MIDAS
%wesite and converts to .mat files as input for the registration process

%imregister() is used to initially align the images using rigid
%registration. The output images are taken as input for the non-rigid
%registration process

%INPUT: filenameI0: the filename of the source image 
%filenameI1: the filename of the target image
%OUTPUT: Moving: matrix of the source image
%Fixed: matrix of the target image

%read .mha image volume files
[BrainI0] = double(mha_read_volume(filenameI1));
[BrainI1] = double(mha_read_volume(filenameI0));
BrainI0 = resize3Dmatrix(size(BrainI1,1),size(BrainI1,2),size(BrainI1,3),BrainI0);

fixed1 = BrainI0;
moving1 = BrainI1;

movingRegistered = zeros(size(BrainI0));
[optimizer,metric] = imregconfig('monomodal');
for i = 1:size(fixed1,3),
    movingRegistered(:,:,i) = imregister(moving1(:,:,i),fixed1(:,:,i),'rigid',optimizer,metric);
end

BrainI1 = movingRegistered;
BrainI0=(BrainI0-min(BrainI0(:)))/(max(BrainI0(:))-min(BrainI0(:)))*11;
BrainI1=(BrainI1-min(BrainI1(:)))/(max(BrainI1(:))-min(BrainI1(:)))*11;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert image to 0-255
Img1=BrainI0*23.1818;
Img2=BrainI1*23.1818;
Img1 = permute(Img1,[1,3,2]);
Img2 = permute(Img2,[1,3,2]);
Moving = Img1; %moving Image
Fixed = Img2;  %fixed Image
end