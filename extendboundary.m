function psi=extendboundary(phi)
% 
psi=phi;
% 
% psi(1,:)=x(1,:);
% psi(end,:)=x(end,:);
% psi(:,1)=y(:,1);
% psi(:,end)=y(:,end);

psi(:,end)=psi(:,end-1);
psi(:,1)=psi(:,2);
psi(1,:)=psi(2,:);
psi(end,:)=psi(end-1,:);


