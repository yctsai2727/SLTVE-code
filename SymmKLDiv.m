function [d]=SymmKLDiv(d1,d2,epsi)
    d1=squeeze(d1);
    d2=squeeze(d2);
    assert(size(d1)==size(d2),"nonconformant arguments in SymmKLDiv");
    if nargin <=2
        epsi = 1e-10;
    end
    d3 = 0.5*(d1+d2);
    d=sqrt(0.5*KLDiv(d1,d3)+0.5*KLDiv(d2,d3));
end

function [d]=KLDiv(d1,d2)
    ldd = log(d1)-log(d2);
    ldd(isnan(ldd)) = 0;
    dldd = d1.*ldd;
    dldd(isnan(dldd)) = 0;
    d=sum(sum(dldd));
end