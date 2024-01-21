function [d]=TotalVariation(d1,d2,l)
    d1=squeeze(d1);
    d2=squeeze(d2);
    assert(size(d1)==size(d2),"nonconformant arguments in TotalVariation");
    if nargin == 2
        l=1;
    end
    d=sum(sum(abs(d1-d2).^l))/2;