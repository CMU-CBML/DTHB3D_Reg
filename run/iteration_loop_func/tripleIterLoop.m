function [pxx, pyy, pzz] = tripleIterLoop(sizeImage, Pixel, Jm, ACP)
%In this function we compute the new positions of the pixel coordinates
%using the spatial transformation function

% INPUT:
% sizeImage: size of the image
% Pixel: the array containing the active element indices and the phi values
% for each pixel
% Jm: the array containing the non-zero splines in each active element
% ACP: the array containing the control points that are set to be active

% OUTPUT:
% pxx, pyy, pzz: the position of the pixel coordinates after spatial transformation

%%
pxx = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
pyy = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));
pzz = zeros(sizeImage(1,1),sizeImage(1,2),sizeImage(1,3));

parfor i = 1:sizeImage(1,3)
    [ temp_pxx, temp_pyy, temp_pzz ] = tripleIterLoopBody(i, sizeImage, Pixel, Jm, ACP  );
    pxx(:,:,i) = temp_pxx;
    pyy(:,:,i) = temp_pyy;
    pzz(:,:,i) = temp_pzz;
end

end

function [ pxx, pyy, pzz ] = tripleIterLoopBody(i, sizeImage, Pixel, Jm, ACP  )

pxx = zeros(sizeImage(1,1),sizeImage(1,2));
pyy = zeros(sizeImage(1,1),sizeImage(1,2));
pzz = zeros(sizeImage(1,1),sizeImage(1,2));

for j = 1:sizeImage(1,2)
    for k = 1:sizeImage(1,1)
        
        % global index of the pixel
        px = sizeImage(1,1)*sizeImage(1,2)*(i-1)+sizeImage(1,1)*(j-1)+k;
        
        % the active element conatining the pixel
        ac_ind = Pixel(px,1).active_cell;
        
        %phi value associated with the pixel coordinate
        supp = Pixel(px,1).phi;
        
        %the non-zero splines over the activ element containing the pixel
        SB = Jm(ac_ind,1).nzsplines;
        
        %control points over the active element
        pts = ACP(SB,1:3);
        
        %compute the new position of the pixel coordinates
        FXX = pts'*supp;
        
        pxx(k,j) = FXX(1,1);
        pyy(k,j) = FXX(2,1);
        pzz(k,j) = FXX(3,1);
    end
end

end

