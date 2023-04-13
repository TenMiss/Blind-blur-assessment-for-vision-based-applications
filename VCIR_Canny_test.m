%%%  This is the main program for blur detection using Canny detector
%
%    Copyright by Wu shiqian 6 March 2006
%    Revised on 20 Oct 2008
clc
clear
close all
%% load images
% [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Select image');
% if isequal([filename,pathname],[0,0])
%     return
% end
% pathAndFilename=strcat(char(pathname),char(filename));
% [pathstr,name,ext,versn] = fileparts(filename);
% I0=imread(pathAndFilename);

I0=imread('C:\My programs\Blur identification\Black & white\parrots_gry.tif');
%I0=imread('C:\Black & white\lake_original.tif');
%I0=imread('J:\Benchmark images\Black & white\houses_gry.tif');
[r c h] = size(I0);
if(h ~= 1)
    I0 = rgb2gray(I0);
end
% figure,imshow(I0);
Tresults =[];
%% Start work
for bsz = 5:2:15

I=I0;
disk = fspecial('disk', bsz);
I = imfilter(I,disk,'circ','conv');

% Pre-defined parameters 
minimumArea = 30;
minimumLineLength =10;
boundary_thresh =30;
minimumAngle =1.5;
LSFLength=21
sigma =1;
thresh =[];
%   Edge detection
[BW0Canny,angleInfo,gmap] = cannyEdgeDetection(I,sigma,thresh);
%figure,imshow(BW0Canny);
%%% Delete noise
BWCanny = bwareaopen(BW0Canny,minimumArea);
[LabelImage,objectNumber]=bwlabel(BWCanny);
imshow(LabelImage)
stats = regionprops(LabelImage,'PixelIdxList','PixelList');
uselines={};
for i =1:objectNumber
    linearInd = stats(i).PixelIdxList;
    localAngles = angleInfo(linearInd)*180/pi;
    possibleLines = lineFind(linearInd,localAngles,minimumLineLength,minimumAngle);    
    if isempty(possibleLines)
        continue
    else
        %%% Determine the line parameters
        for k=1:length(possibleLines)
            BW = repmat(0,r,c);
            rc = possibleLines{k};
            BW(rc)=1;
            BW = bwperim(BW);
            %%%%%%%%%%%%%%%%%%%%%%%% Radon transform
            theta_step=1;
            theta = 0:theta_step:179;
            [R,xp] = radon(BW,theta);
            [rr,cc]=find(R==max(max(R)));
            if length(rr)>1
                disp('This is not a line')
                continue
            else
                line_angle = cc-1; %% Note minus 1 is the real angle
                line_position = rr;
                oneline={[line_angle line_position] rc};
                uselines=cat(1,uselines,oneline);
            end
        end
    end
end
%keyboard
if isempty(uselines)
    disp('Cannot find one suitable line')
else
    %%% Accurately localize the lines
    % [n,m]=size(BWCanny);
    blurParameter1 =zeros(1,size(uselines,1));%blurParameter is a row vector
    blurParameter2 =zeros(1,size(uselines,1));
    I=double(I);
    PSFInfo=[];
    hold on
    for k=1:size(uselines,1)
        objPos = uselines{k,1};
        objPix = uselines{k,2};
        jop=objPos(1);
        iop=objPos(2);
        R_i=jop*theta_step;      %% (R_i is the real angle)
        xp_i=xp(iop);             %%  xp is the distance s to the origin 
        [rr,cc] = ind2sub([r,c],objPix);
        pObj = repmat(0,r,c);
        pObj(objPix)=1;
        %keyboard
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The detected line is divided in 3 cases.
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        if R_i<=1 || (180-R_i)<=1   %% This is a verical line
            drawline(BWCanny,xp_i,R_i)
            pointPositionExtracted =[];
            [x numx]= dataExtraction(cc);            
            if x<=boundary_thresh || (size(I,2)-x)<=boundary_thresh
                disp('Case 1:This vertical line is very near the boarder')
                blurIndex1 = nan;
                blurIndex2 = nan;
            else
                tmp =find(cc==x);
                rowCandidate = rr(tmp);
                for j1=1:length(rowCandidate)
                    extratPoints = x-fix((LSFLength-1)/2):x+fix((LSFLength-1)/2);
                    a =BW0Canny(rowCandidate(j1),extratPoints);
                    
                    [LB,num]=bwlabeln(a);
                    if num ==1
                        pointPositionExtracted =[pointPositionExtracted rowCandidate(j1)];
                        plot(extratPoints,rowCandidate(j1)*ones(1,length(extratPoints)),'r');
                    end
                end
                if length(pointPositionExtracted)<=2
                    disp('Case 1:The vertical line is not available')
                    blurIndex1 = nan;
                    blurIndex2 = nan;                    
                else
                    pointExtracted =I(pointPositionExtracted,x-fix((LSFLength-1)/2):x+fix((LSFLength-1)/2));
                    [LSF,FC] = extractLSF(pointExtracted);
                    if isempty(LSF)
                        disp('Case 1:This line is not suitable which is detected in extractLSF file')
                        blurIndex1 = nan;
                        blurIndex2 = nan;
                    else
                        blurIndex1 = simpleExtractPSF(LSF);
                        blurIndex2 = extractPSF(LSF);
                        drawLine(BWCanny,xp_i,R_i)
                        tempInfo = [R_i,blurIndex1,blurIndex2];
                        PSFInfo = [PSFInfo;tempInfo];
                    end
                end
            end
        elseif abs(R_i-90)<=1     %% This is a horizontal line
            drawline(BWCanny,xp_i,R_i)
            pointPositionExtracted =[];
            [y numy]= dataExtraction(rr);
            if y<=boundary_thresh || (size(I,1)-y)<=boundary_thresh
                disp('Case 2:This horizontal line is very near the boarder')
                blurIndex1 = nan;
                blurIndex2 = nan;
            else
                tmp =find(rr==y);
                columnCandidate = cc(tmp);
                for j1=1:length(columnCandidate)
                    extratPoints = y-fix((LSFLength-1)/2):y+fix((LSFLength-1)/2);
                    a =BW0Canny(extratPoints,columnCandidate(j1));
                    
                    [LB,num]=bwlabeln(a);
                    if num ==1
                        pointPositionExtracted =[pointPositionExtracted columnCandidate(j1)];
                        plot(columnCandidate(j1)*ones(1,length(extratPoints)),extratPoints,'r');
                    end
                end
                if length(pointPositionExtracted)<=2
                    disp('Case 2:The horizonal line is not available')
                    blurIndex1 = nan;
                    blurIndex2 = nan;
                else
                    pointExtracted =I(y-fix((LSFLength-1)/2):y+fix((LSFLength-1)/2),pointPositionExtracted);
                    [LSF,FC] = extractLSF(pointExtracted');
                    if isempty(LSF)
                        disp('Case 2:This line is not suitable which is detected in extractLSF file')
                        blurIndex1 = nan;
                        blurIndex2 = nan;
                    else
                        blurIndex1 = simpleExtractPSF(LSF);
                        blurIndex2 = extractPSF(LSF);
                        drawLine(BWCanny,xp_i,R_i)
                        tempInfo = [R_i,blurIndex1,blurIndex2];
                        PSFInfo = [PSFInfo;tempInfo];
                    end
                end
            end
        else
            %% %%%% For other line
            drawline(BWCanny,xp_i,R_i)
            %%%%%%%%%%  Further find the exact location of the line
            [w,SS]=extractPointsOnLine(pObj,xp_i,R_i);
            lineAxesDetected = w;
            lineCandidate =find(SS~=0);
            pointExtracted =[];
            for j2=1:length(lineCandidate)
                xc =lineAxesDetected(lineCandidate(j2),2);
                yr =lineAxesDetected(lineCandidate(j2),1);
                s = -fix((LSFLength-1)/2):fix((LSFLength-1)/2);
                xxc=xc +s*cos(R_i*pi/180);
                yyr=yr -s*sin(R_i*pi/180);
                roundxc =round(xxc);
                roundyr =round(yyr);
                if min(roundxc) >1 && max(roundxc)< size(I,2) && min(roundyr) >1 & max(roundyr)<size(I,1)
                    a =BW0Canny(roundyr,roundxc);
                    
                    [LB,num]=bwlabeln(a);
                    if num ==1
                        P = interpolateNewValue(I,yyr,xxc);
                        pointExtracted =[pointExtracted;P];
                        plot(roundxc,roundyr,'r');
                    end
                end
            end
            if size(pointExtracted,1)<=2
                disp('Case3:This line is not available')
                blurIndex1 = nan;
                blurIndex2 = nan;
            else
                [LSF,FC] = extractLSF(pointExtracted);
                if isempty(LSF)
                    disp('Case 3:This line is not suitable which is detected in extractLSF file')
                    blurIndex1 = nan;
                    blurIndex2 = nan;
                else
                    blurIndex1 = simpleExtractPSF(LSF);
                    blurIndex2 = extractPSF(LSF);
                    drawLine(BWCanny,xp_i,R_i)
                    tempInfo = [R_i,blurIndex1,blurIndex2];
                    PSFInfo = [PSFInfo;tempInfo];
                end
            end
        end
        PSFInfo
        blurParameter1(k,:)=blurIndex1;
        blurParameter2(k,:)=blurIndex2;
        pause
    end  
end
hold off
if isempty(PSFInfo)
    disp('The blurParameter1 is empty')
else
    LOC = find(PSFInfo(:,2)==min(PSFInfo(:,2)));
    metric_1_minR = PSFInfo(LOC,2)
    UsedLines = length(PSFInfo(:,1));
end
Lresult =[bsz,metric_1_minR];
Tresults =[Tresults;Lresult];
end

Tresults

