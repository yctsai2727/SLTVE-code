function psi=extendboundary(phi)

psi=phi;

% 
%  psi(:,end)=psi(:,end-1);
%  psi(:,1)=psi(:,2);
% psi(1,:)=psi(2,:);
%  psi(end,:)=psi(end-1,:);
 
  psi(1,:)=2*psi(2,:)-psi(3,:);
psi(end,:)=2*psi(end-1,:)-psi(end-2,:);
psi(:,1)=2*psi(:,2)-psi(:,3);
psi(:,end)=2*psi(:,end-1)-psi(:,end-2);
