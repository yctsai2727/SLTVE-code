function [d]=DisplacedSlicedWasserstein(dist1,dist2,x,y,u,k)
    d1=squeeze(dist1{1});
    d2=squeeze(dist2{1});
    assert(size(d1)==size(d2),"nonconformant arguments in SlicedWasserstein");
    [m,n] = size(d1);
    if nargin==5
        k = 5;
    end

    x0 = x(dist1{2},dist1{3});
    y0 = y(dist1{2},dist1{3});
    [d1p,d2p,xp,yp] = Padding(d1,d2,x,y,m,n,x0,y0,u);

    sample = [xp(:) yp(:)];
    P = rand(2,k);
    P = P./vecnorm(P,2);
    Proj = sample*P;

    [Proj_s,P_ind] = sort(Proj,1);
    d1_s = d1p(:)(P_ind);
    d2_s = d2p(:)(P_ind);
    if u==1
        delta = max((Proj_s-[Proj_s(1,:);Proj_s(1:(m+1)*n-1,:)]),([Proj_s(2:end,:);Proj_s((m+1)*n,:)]-Proj_s));
        d = max(sum(abs((d1_s-d2_s).*((m+1)*n:-1:1)).*delta));
    else
        delta = max((Proj_s-[Proj_s(1,:);Proj_s(1:m*(n+1)-1,:)]),([Proj_s(2:end,:);Proj_s(m*(n+1),:)]-Proj_s));
        d = max(sum(abs((d1_s-d2_s).*(m*(n+1):-1:1)).*delta));
    end
end

function [d1p,d2p,xp,yp] = Padding(d1,d2,x,y,m,n,x0,y0,u)
    if u == 1
        d1p = [zeros(m,1),d1(:,:)];
        d2p = [d2(:,:),zeros(m,1)];
        dx = abs(x(2,1)-x(1,1));
        xp = x-x0;
        xp = [xp,xp(m,:)+dx];
        yp = y-y0;
        yp = [yp,yp(m,:)];
    else
        d1p = [zeros(1,n);d1(:,:)];
        d2p = [d2(:,:);zeros(1,n)];
        dy = abs(y(1,2)-y(1,1));
        xp = x-x0;
        xp = [xp;xp(:,n)];
        yp = y-y0;
        yp = [yp;yp(:,n)+dy];
    end
end

    