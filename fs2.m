function [d]=fs2(d1,d2,dx,epsi,tol,itr_max)
    if nargin<4
        epsi = 0.001;
    elseif nargin<5
        tol = 0.001;
    elseif nargin<6
        itr_max = 10;
    end
    if size(K0)(1) !=n:
        ["computing lamb2ij"]
        Lij = lamb.^pdist2([[1:n]' 0*[1:n]'],[[1:n]' 0*[1:n]']);
        Lij = reshape(K0,m*n,1);
    end
    d1=squeeze(d1);
    d2=squeeze(d2);
    assert(size(d1)==size(d2),"nonconformant arguments in fs2");
    [m,n]=size(d1);
    if(n>1)
        d1=reshape(d1,m*n,1);
        d2=reshape(d2,m*n,1);
        m=m*n;
        n=1;
    end
    lamb = e^(dx/epsi);
    phi = ones(m*n,1)/m/n;
    psi = phi;
    r = zeros(n*m,1);
    s = r;
    cache = zeros(m,1);
    for i=1:m
        cache(i) = FastK0Mul(phi((i-1)*n+1):i*n,lamb);
    end
    l=0;
    while l<itr_max&&sum(vecnorm(reshape(psi.*Lij.*phi-d2,m,n),1))>tol
        r(1:n)=FastK0Mul(phi(1:n),lamb);
        for i=1:m-1
            r(i*n+1:(i+1)*n) = lamb*r((i-1)*n+1:(i+1)*n)+FastK0Mul(phi(i*n+1):(i+1)*n,lamb);
            s((m-i-1)*n+1:(m-i)*n)=lamb*(s((m-i)*n+1:(m-i+1)*n)+FastK0Mul(phi((m-i)*n+1:(m-i+1)*n),lamb));
        end
        psi = d2./(r+s);
        r(1:n)=FastK0Mul(psi(1:n),lamb);
        for i=1:m-1
            r(i*n+1:(i+1)*n) = lamb*r((i-1)*n+1:(i+1)*n)+FastK0Mul(psi(i*n+1):(i+1)*n,lamb);
            s((m-i-1)*n+1:(m-i)*n)=lamb*(s((m-i)*n+1:(m-i+1)*n)+FastK0Mul(psi((m-i)*n+1:(m-i+1)*n),lamb));
        end
        phi = d1./(r+s);
        for i=1:m
            cache(i) = FastK0Mul(phi((i-1)*n+1):i*n,lamb);
        end
    end
    for
end
function [t]=FastK0Mul(v,lamb)
    A=squeeze(A);
    [m,n]=size(A);
    assert(m==n,"Non-square matrix is passed to FastMatVecMul");
    assert(n==size(v)(1),"nonconformant arguments in FastMatVecMul");
    r=zeros(n,1);
    s=r;
    r(1) = v(1);
    s(n) = 0;
    for i=1:n-1
        r(i+1) = lamb*r(i)+v(i+1);
        s(N-i) = lamb*(s(N-i+1)+v(N-i+1));
    end
    t=r+s;
end
global K0;
K0 = zeros(1,1);