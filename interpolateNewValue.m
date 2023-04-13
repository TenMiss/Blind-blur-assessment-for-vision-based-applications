function P=interpolateNewValue(I,YR,XC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    When extract LSF,interpolate new value when theta is not zero or 90
%    degree
%
%    YR and XC are computed axes values which are integers
%    Both YR and XC are vectors
%
%    Output P is a row vector (LSF) which contains the pixel values
%
%    Copyright by Wu shiqian 30 Jan 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n =length(YR);
P =zeros(1,n);
for j =1:n
    %%%%%  Interpolate new value
    pxc = XC(j);
    pyr = YR(j);
    Lr = floor(pyr);
    Lc = floor(pxc);
    xcd =pxc-Lc;
    yrd =pyr-Lr;
    p1 =I(Lr,Lc)*(1-yrd)+I(Lr+1,Lc)*yrd;
    p2 =I(Lr,Lc)*(1-xcd)+I(Lr,Lc+1)*xcd;
    p3 =I(Lr,Lc+1)*(1-yrd)+I(Lr+1,Lc+1)*yrd;
    p4 =I(Lr+1,Lc)*(1-xcd)+I(Lr+1,Lc+1)*xcd;
    P(j) = 0.5*(p1*(1-xcd)+p3*xcd +p2*(1-yrd)+p4*yrd);
end