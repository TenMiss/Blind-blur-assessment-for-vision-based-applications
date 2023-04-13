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
tic
[r c h] = size(I0);
if(h ~= 1)
    I0 = rgb2gray(I0);
end
% figure,imshow(I0);
I=I0;
disk = fspecial('disk', 3);
I = imfilter(I,disk,'circ','conv');
%figure,imshow(I);

%%%   Edge detection
[BW0Sobel,angleInfo] = SobelEdgeDetection(I);
BW0Canny = edge(I,'canny');
%%% Delete noise
BWSobel = bwareaopen(BW0Sobel,5);
%BWCanny = bwareaopen(BW0Canny,5);
figure,imshow(BWSobel);
%figure,imshow(BWCanny);
[LabelImage,objectNumber]=bwlabel(BWSobel);
%[LabelImage,objectNumber]=bwlabel(BWCanny);
stats = regionprops(LabelImage,'all');
uselines={};
for i =1:objectNumber
    linearInd = stats(i).PixelIdxList;
    localGrad = angleInfo(linearInd)*180/pi;
    medianGrad = median(localGrad);
    numObj=length(localGrad);
    possibleLines={};
    while numObj>10
        [x,localGrad,num]=lineFind(localGrad,medianGrad,5);
        if isempty(x)
            delt = localGrad-medianGrad;
            num = find(delt==min(delt));
            localGrad(num)=[];
            linearInd(num)=[];        
        elseif length(x)>=10
            pos = linearInd(num);
            possibleLines = cat(1,possibleLines,{pos});
            linearInd(num)=[];
        else
            linearInd(num)=[];        
        end
        numObj=length(localGrad);
        medianGrad = median(localGrad);
    end
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
    % [n,m]=size(BWSobel);
    blurParameter1 =zeros(1,size(uselines,1)); % i.e., blurParameter is a 1*row vector
    blurParameter2 =zeros(1,size(uselines,1));
    I=double(I);
    boundary_thresh =35;
    LSFLength = 31;
    LSFInfo={};
    angle_error = 20
    hold on
    for k=1:size(uselines,1)
        objPos = uselines{k,1};
        objPix = uselines{k,2};
        j=objPos(1);
        i=objPos(2);
        R_i=j*theta_step;      %% (R_i is the real angle)
        angle = R_i;
        xp_i=xp(i);             %%  xp is the distance s to the origin 
        [rr,cc] = ind2sub([r,c],objPix);
        pObj = repmat(0,r,c);
        pObj(objPix)=1;
                
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The detected line is divided in 3 cases.
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        if R_i<=1 | (180-R_i)<=1   %% This is a verical line
            %drawLine(BWSobel,xp_i,R_i)
            pointPositionExtracted =[];
            x=floor(size(I,2)/2 + 0.5) + xp_i;
            if x<=boundary_thresh | (size(BWSobel,2)-x)<=boundary_thresh
                disp('Case 1:This vertical line is very near the boarder')
                blurIndex1 = nan;
                blurIndex2 = nan;
            else
                rowCandidate =find(cc==x);
                for j1=1:length(rowCandidate)
                    extratPoints = x-(LSFLength-1)/2:x+(LSFLength-1)/2;
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
                    pointExtracted =I(pointPositionExtracted,x-(LSFLength-1)/2:x+(LSFLength-1)/2);
                    [LSF,FC] = extractLSF(pointExtracted);
                    if isempty(LSF)
                        disp('Case 1:This line is not suitable which is detected in extractLSF file')
                        blurIndex1 = nan;
                        blurIndex2 = nan;
                    else
                        blurIndex1 = simpleExtractPSF(LSF);
                        blurIndex2 = extractPSF(LSF);
                        drawLine(BWSobel,xp_i,R_i)
                        tempInfo = { blurIndex1,blurIndex2,LSF};
                        LSFInfo = cat(1,LSFInfo,tempInfo);
                    end
                end
            end
        elseif abs(R_i-90)<=1     %% This is a horizontal line
            %drawLine(BWSobel,xp_i,R_i)
            pointPositionExtracted =[];
            y =floor(size(I,1)/2 + 0.5)-xp_i;
            if y<=boundary_thresh | (size(BWSobel,1)-y)<=boundary_thresh
                disp('Case 2:This horizontal line is very near the boarder')
                blurIndex1 = nan;
                blurIndex2 = nan;
            else
                columnCandidate =find(rr==y);
                for j1=1:length(columnCandidate)
                    extratPoints = y-(LSFLength-1)/2:y+(LSFLength-1)/2;
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
                    pointExtracted =I(y-(LSFLength-1)/2:y+(LSFLength-1)/2,pointPositionExtracted);
                    [LSF,FC] = extractLSF(pointExtracted');
                    if isempty(LSF)
                        disp('Case 2:This line is not suitable which is detected in extractLSF file')
                        blurIndex1 = nan;
                        blurIndex2 = nan;
                    else
                        blurIndex1 = simpleExtractPSF(LSF);
                        blurIndex2 = extractPSF(LSF);
                        drawLine(BWSobel,xp_i,R_i)
                        tempInfo = { blurIndex1,blurIndex2,LSF};
                        LSFInfo = cat(1,LSFInfo,tempInfo);
                    end
                end
            end
        else
            %% %%%% For other line
            %drawLine(BWSobel,xp_i,R_i)
            %%%%%%%%%%  Further find the exact location of the line
            [w,SS]=extractPointsOnLine(pObj,xp_i,R_i);
            lineAxesDetected = w;
            lineCandidate =find(SS~=0);
            pointExtracted =[];
            for j2=1:length(lineCandidate)
                xc =lineAxesDetected(lineCandidate(j2),2);
                yr =lineAxesDetected(lineCandidate(j2),1);
                s = -(LSFLength-1)/2:(LSFLength-1)/2;
                xxc=xc +s*cos((R_i+angle_error)*pi/180);
                yyr=yr -s*sin((R_i+angle_error)*pi/180);
                roundxc =round(xxc);
                roundyr =round(yyr);
                if min(roundxc) >1 & max(roundxc)< size(I,2) & min(roundyr) >1 & max(roundyr)<size(I,1)
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
                    drawLine(BWSobel,xp_i,R_i)
                    tempInfo = { blurIndex1,blurIndex2,LSF};
                    LSFInfo = cat(1,LSFInfo,tempInfo);
                end
            end
        end
        blurParameter1(k)=blurIndex1;
        blurParameter2(k)=blurIndex2;
    end  
end
hold off
if isempty(blurParameter1)
    disp('The blurParameter1 is empty')
else
    non_nan = ~isnan(blurParameter1);
    n = find(non_nan~=0);
    parameters11=blurParameter1(n);
    if isempty(n)
        disp('There are no suitable parameter in blur parameter')
    else
         metric_1_minR = min(blurParameter1(n))
         metric_1_medianR = median(blurParameter1(n));
    end
end
time=toc





