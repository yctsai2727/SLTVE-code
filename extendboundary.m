function psi=extendboundary(phi)
    psi=phi;
    psi(:,end)=psi(:,end-1);
    psi(:,1)=psi(:,2);
    psi(1,:)=psi(2,:);
    psi(end,:)=psi(end-1,:);


