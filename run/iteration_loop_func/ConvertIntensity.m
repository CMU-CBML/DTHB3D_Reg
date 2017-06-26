function [ Ilabel ] = ConvertIntensity( I0 )
%source code taken from the website https://github.com/stellaccl/cdmffd-image-registration
%Paper: Chan, Chiu Ling, et al. "Two and Three Dimensional Image Registration Based on B-Spline Composition
% and Level Sets." Communications in Computational Physics 21.2 (2017): 600-622.
%ConverIntensity convert intensity 0-255 to label 0-11
% input I0= 3 dimensional image
% Output Ilabel = Image with label 0-11

nx=size(I0,2);
ny=size(I0,1);
nz=size(I0,3);

Ilabel=zeros(size(I0));

for i=1:nx
    for j=1:ny
        for k=1:nz
            if I0(j,i,k)>=0 && I0(j,i,k)<=21.25
                Ilabel(j,i,k)=0;
            end
            if I0(j,i,k)>21.25 && I0(j,i,k)<=42.5
                Ilabel(j,i,k)=1;
            end
            if I0(j,i,k)>42.5 && I0(j,i,k)<=63.75
                Ilabel(j,i,k)=2;
            end
            if I0(j,i,k)>63.75 && I0(j,i,k)<=85
                Ilabel(j,i,k)=3;
            end
            if I0(j,i,k)>85 && I0(j,i,k)<=106.25
                Ilabel(j,i,k)=4;
            end
            if I0(j,i,k)>106.25 && I0(j,i,k)<=127.5
                Ilabel(j,i,k)=5;
            end
            if I0(j,i,k)>127.5 && I0(j,i,k)<=148.75
                Ilabel(j,i,k)=6;
            end
            if I0(j,i,k)>148.75 && I0(j,i,k)<=170
                Ilabel(j,i,k)=7;
            end
            if I0(j,i,k)>170 && I0(j,i,k)<=191.25
                Ilabel(j,i,k)=8;
            end
            if I0(j,i,k)>191.25 && I0(j,i,k)<=212.5
                Ilabel(j,i,k)=9;
            end
            if I0(j,i,k)>212.5 && I0(j,i,k)<=233.75
                Ilabel(j,i,k)=10;
            end
            if I0(j,i,k)>233.75 && I0(j,i,k)<=255
                Ilabel(j,i,k)=11;
            end
        end
    end
end
end

