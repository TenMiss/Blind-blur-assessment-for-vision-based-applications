clear all
clc
R=20;
x=-R:R;
y=2*sqrt(R*R-x.*x)/(pi*R*R);
R=20;
x1=-R:R;
y1=2*sqrt(R*R-x1.*x1)/(pi*R*R);

delt =15
d=zeros(1,delt)
xx = -20:35;
yy = [y d];
yy1 = [d y1];
figure,plot(xx,yy,'r-',xx,yy1,'k-')
axis([-30 40 0 0.04])
xlabel('Position')
ylabel('LSF')

figure,plot(xx,yy,'r--',xx,yy1,'r--')
hold on
plot(xx,yy+yy1,'k-')
axis([-30 40 0 0.06])
xlabel('Position')
ylabel('LSF')
hold off
