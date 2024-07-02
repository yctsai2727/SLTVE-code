function [d]=DisplacedTotalVariation(dist1,dist2,dirc,l)
    d1=squeeze(dist1{1});
    d2=squeeze(dist2{1});
    [m,n]=size(d1);
    assert(size(d1)==size(d2),"nonconformant arguments in TotalVariation");
    if nargin == 3
        l=1;
    end
    [d1p,d2p] = Padding(d1,d2,m,n,dirc);
    d=sum(sum(abs(d1p-d2p)));
end

function [d1p,d2p] = Padding(d1,d2,m,n,u)
    if u == 1
        d1p = [zeros(1,n);d1(:,:)];
        d2p = [d2(:,:);zeros(1,n)];
    else
        d1p = [zeros(m,1),d1(:,:)];
        d2p = [d2(:,:),zeros(m,1)];
    end
end