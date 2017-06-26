function [CellGrad, meanGrad] = computeDiffGradImage(cell_co,I1,I2, pixX, pixY, pixZ)

% In this function the absolute value of the image difference between the
% evolving image and target image.

% INPUT:
% cell_co: the coordinates of the centroid of the active elements
% I1: source image
% I2: target image
% pixX, pixY, pixZ: the coordinates of the pixels of the image

% OUTPUT:
% CellGrad: the 
% meanGrad: 

%#codegen
Idiff = abs(I1-I2);
%[DDI1X,DDI1Y,DDI1Z] = gradient((I1-I2));
%[DDDI1X,DDDI1Y,DDDI1Z] = gradient((Img1-Img2));
%[DDI1X,DDI1Y,DDI1Z] = gradient((I1));
%Idiff = sqrt((DDI1X).^2 + (DDI1Y).^2 + (DDI1Z).^2);
%Idiff_in = sqrt((DDDI1X).^2 + (DDDI1Y).^2 + (DDDI1Z).^2);
Cell_grad = interp3(pixY,pixX,pixZ,Idiff,cell_co(:,2),cell_co(:,1),cell_co(:,3));

meanGrad= mean(mean(mean(Idiff,1),2),3);
CellGrad = Cell_grad;

end