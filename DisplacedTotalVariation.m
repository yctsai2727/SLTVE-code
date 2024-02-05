function [d]=DisplacedTotalVariation(d1,d2,dirc,l) #dirc=1: horizontal shifting, dirc=2: vertical shifting
    d1=squeeze(d1);
    d2=squeeze(d2);
    [m,n]=size(d1);
    assert(size(d1)==size(d2),"nonconformant arguments in DisplacedTotalVariation");
    if nargin == 3
        l=1;
    end
    if dirc==1
        d=TotalVariation(d1(:,1:n-1),d2(:,2:n),l);
    elseif
        d=TotalVariation(d1(1:m-1,:),d2(2:m,:),l);
    end