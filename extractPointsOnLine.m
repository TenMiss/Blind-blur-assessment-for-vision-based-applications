function [w,SBW]=extractPointsOnLine(BW,s,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Extract points from BW for a specific line
%
%    Input parameters s and beta are obtained from radon transform of BW
%
%    Output w and SBW are axes of the detected line and the corresponding
%    valeus of BW
%    w is m by 2 matrix where first column is the axes of row (y) and second column
%    is the axes of column
%    Note: x is in ascending order but y may not.
%    SBW is a column vector containing m elements which shows whether the line is on 
%    where [n,m]=size(BW)
%
%    Copyright by Wu shiqian 30 Jan 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,m]=size(BW);
x_origin=m/2 +s*cos(beta*pi/180);
y_origin=n/2 -s*sin(beta*pi/180);
xv=1:m;
yv=(y_origin-(xv-x_origin)*tan((beta-90)*pi/180));
w=[yv' xv'];
w=round(w);
SBW=zeros(m,1);
for j=1:m
    if w(j,1)<1 | w(j,1)> n
        SBW(j)=0;
    else
        SBW(j)= BW(w(j,1),w(j,2));
    end
end