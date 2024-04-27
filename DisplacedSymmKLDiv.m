function [d]=DisplacedSymmKLDiv(d1,d2,u,epsi)
    d1=squeeze(d1);
    d2=squeeze(d2);
    [m,n] = size(d1);
    if u == 1:
        d1 = d1(2:m,:);
        d2 = d2(1:m-1,:);
    else:
        d1 = d1(:,2:n);
        d2 = d2(:,1:n-1);
    end
    d3=0.5*(d1+d2);
    assert(size(d1)==size(d2),"nonconformant arguments in DisplacedSymmKLDiv");
    d=SymmKLDiv(d1,d2,epsi);