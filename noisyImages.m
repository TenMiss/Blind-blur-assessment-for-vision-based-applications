%%%  This is the main program for blur detection
%
%    Copyright by Wu shiqian 6 March 2006
clc
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
if ~isa(I0, 'double')
  I0 = im2double(I0);
  classChanged = 1;
end
% figure,imshow(I0);
I=I0;
disk = fspecial('disk', 8);
I = imfilter(I,disk,'circ','conv');    
figure,imshow(I);

IV=I(:);
Imean=mean(IV);
I1= IV-Imean;
ICov = sum(I1.*I1)/(r*c)
SNR=40
noiseCov=ICov./(10.^(SNR/10))
It = imnoise(I,'gaussian',0,noiseCov);
figure,imshow(It);

nI=I-It;
VI=nI(:);
NVan = sum((VI-mean(VI)).*(VI-mean(VI)))/(r*c)



