function psi=fixboundary(phi,x)

psi=phi;

psi(1,:)=x(1,:);
psi(end,:)=x(end,:);
psi(:,1)=x(:,1);
psi(:,end)=x(:,end);

end