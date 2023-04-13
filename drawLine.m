function drawLine(BW,s,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Draw a line on the detected line to confirm whether it is correct
%
%    Parameters s and beta are obtained from radon transform
%
%    Copyright by Wu shiqian 30 Jan 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,m]=size(BW);
x_origin=m/2 +s*cos(beta*pi/180);
y_origin=n/2 -s*sin(beta*pi/180);
x1=1;
xe=m;
y1=(y_origin-(x1-x_origin)*tan(((beta)-90)*pi/180));
ye=(y_origin-(xe-x_origin)*tan(((beta)-90)*pi/180));
xv=[x1 xe];
yv=[y1 ye];
line(xv,yv)
plot(x_origin, y_origin,'.r')