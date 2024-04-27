function [d]=DisplacedTotalVariation(d1,d2,dirc,l)
    d1=squeeze(d1);
    d2=squeeze(d2);
    [m,n]=size(d1);
    assert(size(d1)==size(d2),"nonconformant arguments in TotalVariation");
    if nargin == 3
        l=1;
    end
    if dirc==2
        d=sum(sum(abs(d1(:,1:n-1)-d2(:,2:n)).^l));
    elseif
        d=sum(sum(abs(d1(1:m-1,:)-d2(2:m,:)).^l));
    end