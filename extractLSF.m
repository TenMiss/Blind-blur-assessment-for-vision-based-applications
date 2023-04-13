function [FLSF,FC] = extractLSF(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Extract LSF from the area of the original image which contains a line
%    Note:
%%   Input A: Each row represents a sample     
%
%    Copyright by Wu shiqian 30 Jan 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
centerPosition = floor(size(A,2)/2 +0.5);
rowmax =max(A,[],2);
maxMatrix =rowmax * ones(1,size(A,2));
A = A ./maxMatrix;
LSF=abs(diff(A,1,2 ));
LSF(:,end+1)=0;
aLSF=[];
FC=[];
for i=1:size(LSF,1)
    a=LSF(i,:);
    [amax,CC]=max(a);
    if abs(CC-centerPosition)<=3
        LSF1 =circshift(a,[0,centerPosition-CC]);
        aLSF = [aLSF;LSF1];
        [amax2,CC2]=max(LSF1);
        FC=[FC CC2];
    end
end
if size(aLSF,1)<=2
    %disp('This is not a suitable line which is detected in extractLSF file')
    FLSF =[];
else
    FLSF=median(aLSF);
end


