function [d]=SymmKLDiv(d1,d2,epsi)
    d1=squeeze(d1);
    d2=squeeze(d2);
    d3=0.5*(d1+d2);
    assert(size(d1)==size(d2),"nonconformant arguments in SymmKLDiv");
    if nargin <=2
        epsi = 1e-10;
    end
    d=sqrt(KLDiv(d1,d3,epsi)+KLDiv(d2,d3,epsi));

function [d]=KLDiv(d1,d2,epsi)
    d2=d2.*(d2>0)+epsi*(d2<=0);
    d=sum(sum(d1.*(log(d1)-log(d2))));