function [ pxx, pyy, pzz ] = tripleIterLoopBody(i, sizeImage, Pixel, Jm, ACP  )
pxx = zeros(sizeImage(1,1),sizeImage(1,2));
pyy = zeros(sizeImage(1,1),sizeImage(1,2));
pzz = zeros(sizeImage(1,1),sizeImage(1,2));
for j = 1:sizeImage(1,2)
    for k = 1:sizeImage(1,1)
        px = sizeImage(1,1)*sizeImage(1,2)*(i-1)+sizeImage(1,1)*(j-1)+k;
        ac_ind = Pixel(px,1).active_cell;
        supp = Pixel(px,1).phi;
        SB = Jm(ac_ind,1).nzsplines;
        
        pts = ACP(SB,1:3);
        FXX = pts'*supp;
        pxx(k,j) = FXX(1,1);
        pyy(k,j) = FXX(2,1);
        pzz(k,j) = FXX(3,1);
    end
end
end

