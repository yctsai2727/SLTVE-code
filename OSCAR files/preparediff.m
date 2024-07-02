function [xm,xp,ym,yp]=preparediff(phi)

[m,n]=size(phi);
xm=[phi(1,:);phi(1:m-1,:)];
xp=[phi(2:m,:);phi(m,:)];  
ym=[phi(:,1),phi(:,1:n-1)];
yp=[phi(:,2:n),phi(:,n)]; 
