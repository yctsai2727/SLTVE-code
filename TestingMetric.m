function [d]=TestingMetric(d1,d2,x,y,dx,dy)
    d1=squeeze(d1);
    d2=squeeze(d2);
    assert(size(d1)==size(d2),"nonconformant arguments in TestingMetric");
    [m,n]=size(d1);
    d=sqrt((sum(((d1-d2).*x)(:))*dx*dy)^2+(sum(((d1-d2).*y)(:))*dx*dy)^2);
    