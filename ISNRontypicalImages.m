%%%  This is the main program for blur detection
%    Copyright by Wu shiqian 6 March 2006
clear
close all

[filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Select image');
if isequal([filename,pathname],[0,0])
    return
end
pathAndFilename=strcat(char(pathname),char(filename));
[pathstr,name,ext,versn] = fileparts(filename);
I0=imread(pathAndFilename);
[r c h] = size(I0);
if(h ~= 1)
    I0 = rgb2gray(I0);
end
% figure,imshow(I0);
I=I0;
par =13
disk = fspecial('disk', par-1);
I = imfilter(I,disk,'circ','conv');
%restoration
PSF1 = fspecial('disk',par);
AA1 = edgetaper(I,PSF1);
A1 = deconvlucy(AA1,PSF1);
I0=double(I0);
I=double(I);
A1=double(A1);
sum1= sum(sum((I0-I).^2));
sum2= sum(sum((I0-A1).^2));
ISNR1=10*log10(sum1/sum2)
