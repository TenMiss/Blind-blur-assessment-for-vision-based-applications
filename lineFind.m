function lineForm = lineFind(linearInd,localAngles,minimumLineLength,minimumAngle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      
%               Find lines
%
%  Parameters:
%
%  inputs:     'linearInd' is the position vector.
%              'localAngles' is the angle vector 
%
%  outputs:   'lineForm' is the cell containing the position data which form lines
%
%  Copright by Wu Shiqian  27 Oct 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lineForm = {};
if length(localAngles)<=minimumLineLength
    return
else
    sig=1;
end
while sig
    a = localAngles*ones(1,length(linearInd));
    delt = abs(a-a');
    delt = delt<=minimumAngle;
    TLength = sum(delt);
    [DataLength,dataN]=max(TLength);
    if DataLength >= minimumLineLength
        D = delt(:,dataN);
        n = find(D==1);
        aLine = linearInd(n);
        linearInd(n)=[];
        localAngles(n)=[];
        lineForm = cat(1,lineForm,{aLine});
        if length(linearInd)> minimumLineLength
            sig = 1;
        else
            sig=0;
        end
    else
        sig =0;
    end
end
    
end
