function blurIndex = extractPSF(LSF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Extract PSF from the LSF
%
%    Copyright by Wu shiqian 30 Jan 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(LSF);
centroid=0;
var = 0;
for i=1:n
    centroid = centroid+i*LSF(i);
end
centroid = centroid/sum(LSF);
for i=1:n
    var = var + (i-centroid)^2*LSF(i);   
end
blurIndex=2*sqrt(var);
if blurIndex <1
    blurIndex=nan;
end