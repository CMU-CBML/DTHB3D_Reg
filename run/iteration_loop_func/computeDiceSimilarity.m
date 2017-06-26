function computeDiceSimilarity(fid_out,iterct,I1,I2)
% source code taken from the website https://github.com/stellaccl/cdmffd-image-registration
% Paper: Chan, Chiu Ling, et al. "Two and Three Dimensional Image Registration 
% Based on B-Spline Composition and Level Sets." Communications in Computational Physics 21.2 (2017): 600-622.
% this function computes the dice similarity for the images taken from
% Brainweb website

[ I0label ] = ConvertIntensity( I1 );
[ I1label ] = ConvertIntensity( I2);

LabelSize=12;

nx=size(I1,2);
ny=size(I1,1);
nz=size(I1,3);

II0=I0label(:)';
II1=I1label(:)';

[ Dice ] = DiceSimilarity( nx, ny, nz, II0, II1 ,LabelSize);

fprintf(fid_out, '%d, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f\r\n',iterct, Dice(1), Dice(2), Dice(3), Dice(4), Dice(5), Dice(6),Dice(7), Dice(8), Dice(9), Dice(10), Dice(11), Dice(12));

end