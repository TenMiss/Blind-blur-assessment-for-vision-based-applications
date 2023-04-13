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
% figure,imshow(I0);
I=I0;
for i=1:4
    disk = fspecial('disk', 2*i+3);
    I = imfilter(I,disk,'circ','conv');
    figure,imshow(I);
end




