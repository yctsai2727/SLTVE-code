function [d]=DisplacedSymmKLDiv(dist1,dist2,u,epsi)
    d1=squeeze(dist1{1});
    d2=squeeze(dist2{1});
    assert(size(d1)==size(d2),"nonconformant arguments in DisplacedSymmKLDiv");
    [m,n] = size(d1);
    [d1p,d2p] = Padding(d1,d2,m,n,u);
    if nargin == 3
        d=SymmKLDiv(d1p,d2p);
    else
        d=SymmKLDiv(d1p,d2p,epsi);
    end
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