function [phi,psi,x,y,dx,dy]=initial(m,n,xmin,xmax,ymin,ymax)

dx=(xmax-xmin)/(m-1);
dy=(ymax-ymin)/(n-1);

[x,y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
x=x';
y=y';

phi=x;  
psi=y;
