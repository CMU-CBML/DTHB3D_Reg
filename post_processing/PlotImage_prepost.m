function PlotImage_prepost(iterct,Img1,Img2,Iplot)
% In this function the image is displayed and stored in .png format. 

% INPUT:
% iterct: iteration count
% Img1: source or moving image
% Img2: target image
% Iplot: registered image to be displayed

%%

if(iterct==1)
    figure('Position',[100,100,1300,900])
end

nx = size(Iplot,1);
ny = size(Iplot,2);
nz = size(Iplot,3);

position_slice_x1=round(nx/2);
position_slice_y1=64;
position_slice_z1=round(nz/2);

position_slice_x2=round(nx/2);
position_slice_y2=79;
position_slice_z2=round(nz/2);
% Source image at the start of registration
% slice in y direction
for i=1:nz
    I0x(:,i)=Img1(:,position_slice_y1,i);
end

I0x=I0x';
I0x=flipud(I0x);
subtightplot(3,3,1,[0.05,0.000005])
imagesc(I0x)
title(['I0 at x=', num2str(position_slice_x1),'  iter = ',num2str(iterct)]);
axis image
colormap('gray')
drawnow

%slice at x position
for i=1:nz
    I0y(:,i)=Img1(position_slice_x1,:,i);
end

I0yy=I0y';
I0yy=flipud(I0yy);
subtightplot(3,3,4,[0.05,0.000005])
imagesc(I0yy)
title(['I0 at y= ', num2str(position_slice_y1),'  iter = ',num2str(iterct)])
axis image
colormap('gray')
drawnow

%slice at z position
subtightplot(3,3,7,[0.05,0.000005])
imagesc(flipud(Img1(:,:,position_slice_z1)'))
title(['I0 at z= ', num2str(position_slice_z1),'  iter = ',num2str(iterct)])
axis image
colormap('gray')
drawnow

% Registered image
% slice at x position
for i=1:nz
    I0x(:,i)=Iplot(:,position_slice_y2,i);
end

I0x=I0x';
I0x=flipud(I0x);

subtightplot(3,3,2,[0.05,0.000005])
imagesc(I0x)
title(['I1 at x=', num2str(position_slice_x2),'  iter = ',num2str(iterct)]);
axis image
colormap('gray')
drawnow

%slice at y position
for i=1:nz
    I0y(:,i)=Iplot(position_slice_x2,:,i);
end

I0yy=I0y';
I0yy=flipud(I0yy);
subtightplot(3,3,5,[0.05,0.000005])
imagesc(I0yy)

title(['I1 at y= ', num2str(position_slice_y2),'  iter = ',num2str(iterct)])
axis image
colormap('gray')
drawnow

%slice at z position
subtightplot(3,3,8,[0.05,0.000005])
imagesc(flipud(Iplot(:,:,position_slice_z2)'))
title(['I1 at z= ', num2str(position_slice_z2),'  iter = ',num2str(iterct)])
axis image
colormap('gray')
drawnow

% Target image
% slice at y position
for i=1:nz
    I1x(:,i)=Img2(:,position_slice_y2,i);
end

I1x=I1x';
I1x=flipud(I1x);
subtightplot(3,3,3,[0.05,0.000005])
imagesc(I1x)
title(['I2 at x= ', num2str(position_slice_x2),'  iter = ',num2str(iterct)])
axis image
colormap('gray')
drawnow

%slice at x position
for i=1:nz
    I1y(:,i)=Img2(position_slice_x2,:,i);
end

I1yy=I1y';
I1yy=flipud(I1yy);
subtightplot(3,3,6,[0.05,0.000005])
imagesc(I1yy)
title(['I2 at y= ', num2str(position_slice_y2),'  iter = ',num2str(iterct)])
axis image
colormap('gray')
drawnow

%slice at z position
subtightplot(3,3,9,[0.05,0.000005])
imagesc(flipud(Img2(:,:,position_slice_z2)'))
title(['I2 at z= ', num2str(position_slice_z2),'  iter = ',num2str(iterct)])
axis image
colormap('gray')
drawnow

% save image
filename_img = sprintf('../PostProcessing/evolve_image%d.png',iterct);
saveas(gcf,filename_img)

end

