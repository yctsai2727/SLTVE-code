clear all;
close all;
format compact
format long

load u_2014.mat   
load v_2014.mat  % These two datasets contain the velocity data in the whole year 2014.

starttime=cputime;

[u_field,v_field]=DataPreprocessing(u_2014,v_2014);  % This function is
                                                       % used to extract
                                                       % and convert
                                                       % data we are interested in.
% load u_field.mat;
% load v_field.mat;


%%%%%%%%Discretization%%%%%%%%%%%%%
    [m,n]=size(u_field(:,:,1))
    xmin=0;
    xmax=50;
    ymin=0;
    ymax=25; % each length unit corresponds to one longitude or latitude degree
    finalt=50; % each temporal unit corresponds to one day
%%%%%Discretization%%%%%%%%%%%%%

%%%%%%%Initialization%%%%%%%%%%%
% cfl=0.75;
% [phi,psi,x,y,dx,dy]=initial(m,n,xmin,xmax,ymin,ymax);
% dt=0.25;
% t=0;
% i=0;

% %%%%%%%Initialization%%%%%%%%%%%

% while (t<finalt)
   
%    i=i+1;
	
%   [t,dt,t+dt]   
   
%    %%%%%%%%%%%%Computing the forward flow map \Phi_n^{n+1}%%%%%%%%%%%%%%%%%%%%%
%    u=u_field(:,:,i);v=v_field(:,:,i);  
%    u1=u_field(:,:,i+1);v1=v_field(:,:,i+1);  
%    [xm1, xp1, ym1, yp1] = preparediff(u1);
%    [up, um, vp, vm] = WENO2(u1, xm1, xp1, ym1, yp1, dx, dy);
%    phi=x+0.5*dt*(u+u1)+0.5*dt^2*(u.*up+v.*vp);
%    [xm1, xp1, ym1, yp1] = preparediff(v1);
%    [up, um, vp, vm] = WENO2(v1, xm1, xp1, ym1, yp1, dx, dy);
%    psi=y+0.5*dt*(v+v1)+0.5*dt^2*(u.*up+v.*vp); 
%    phi=fixboundary(phi,x);
%    psi=fixboundary(psi,y);
%    %%%%%%%%%%%%Computing the forward flow map \Phi_n^{n+1}%%%%%%%%%%%%%%%%%%%%%
   
%    %%%%%%%%%%%Interpolate to obtain forward flow map \Phi_0^{n+1}%%%%%%%%%%%   
%    if(t>0)
%     phi=interp2(x',y',phi',phitemp',psitemp','cubic')';
%     psi=interp2(x',y',psi',phitemp',psitemp','cubic')';
%     phi=min(max(phi,xmin),xmax);
%     psi=min(max(psi,ymin),ymax);
%    end   
%    phitemp=phi;
%    psitemp=psi;   
%    %%%%%%%%%%%%%%%%%%Interpolate to obtain forward flow map \Phi_0^{n+1}%%%%%%%%%%%%%%%
%   t=t+dt;
% end

% %%%%%%%%%Compute the FTLE%%%%%%%%%%%%%%%%%
%    [phi12,phi11]=gradient(phi,dy,dx);
%    [phi22,phi21]=gradient(psi,dy,dx);
%    a=phi11.^2+phi21.^2;
%    b=phi11.*phi12+phi21.*phi22;
%    c=phi12.^2+phi22.^2;
%    lambda=(a+c+sqrt((a+c).^2-4*(a.*c-b.^2)))/2; 
%    FTLE=log(max(1e-30,lambda))/2/abs(finalt);
%    imagesc(x(:,1),y(1,:),FTLE',[-0.05,0.1])
%    axis equal
%    axis xy
%    axis([0 50 0 25])
%    colorbar
% %%%%%%%%%Compute the FTLE%%%%%%%%%%%%%%%%%

% runtime=cputime-starttime;