function [Moving,Fixed] = read3DImage_mrT1(filenameI0,filenameI1)

%This function reads the .mha files files of brain MR data taken from RIRE data base (http://www.insight-journal.org/rire/download_data.php) 
% to .mat files as input for the registration process 
%INPUT: filenameI0: the filename of the source image 
%filenameI1: the filename of the target image
%OUTPUT: Moving: matrix of the source image
%Fixed: matrix of the target image
% Image is padded with zeros on the boundaries to make the interpolation
% exact


%read raw binary files
[BrainI0] = double(mha_read_volume(filenameI0));
[BrainI1] = double(mha_read_volume(filenameI1));
BrainI1(:,:,52) = zeros(size(BrainI1,1),size(BrainI1,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert image to 0-255
BrainI0=(BrainI0-min(BrainI0(:)))/(max(BrainI0(:))-min(BrainI0(:)))*255;
BrainI1=(BrainI1-min(BrainI1(:)))/(max(BrainI1(:))-min(BrainI1(:)))*255;

Moving = permute(BrainI0,[1,3,2]); %moving Image
Fixed = permute(BrainI1,[1,3,2]);  %fixed Image
end