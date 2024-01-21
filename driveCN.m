#! octave -qf
clear all
close all
format compact
format long

starttime = cputime;

arglist=argv();
id = str2num(argv{1})

%parameter
m = 129;
n = (m - 1) / 2;
D0 = 0.0001;
method = @TotalVariation;

%domain info
xmin = 0;
xmax = 2;
ymin = 0;
ymax = 1;
finalt = 10;
K=100;

%parallel frame parameter
m_p = 33;
n_p = 32;
id_x=mod(id,4);
id_y=floor(id/4);

%rank approximation
r=15;
decoder = @(M) LowRankDecoder(M,r,m,n);

%initialization
[x, y, dx, dy] = initial(m, n, xmin, xmax, ymin, ymax);
dt0 = 0.05;
dt = 0.05;
dtK = (finalt-dt0)/(K-1);
miu = D0 * dt / dx^2;
M1 = left_mat(miu, m, n);
M2 = right_mat(miu, m, n);
traj = zeros(m_p,n_p,K,r,m+n+1);
cache = zeros(m+1,n+1,2); #1->up, 2->left
flag = zeros(m,n);
sol = zeros(m_p,n_p);

%main loop
for i = 1:m_p
    for j = 1:n_p
        [i, j]
        if id_x*m_p+i>m || id_y*n_p+j>n
            break
        end
        x0 = x(id_x*m_p+i, 1);
        y0 = y(1, id_y*n_p+j);
        t = dt0;
        down=0;
        right=0;
        if flag(i,j)==0
            pdf = 1 / (4 * pi * D0 * dt0) * exp(-1 / (4 * D0 * dt0) * ((x - x0).^2 + (y - y0).^2));
            traj(i,j,:,:,:)=Solver(x,y,x0,y0,t,finalt,dx,dy,dt,dtK,pdf,m,n,K,r,M1,M2);
            %size(traj(i,j,:,:,:))
            flag(i,j)=1;
        end
        if i>1 && flag(i-1,j)==0
            Lpdf = (1 / (4 * pi * D0 * dt0) * exp(-1 / (4 * D0 * dt0) * ((x - x0+dx).^2 + (y - y0).^2)));
            traj(i-1,j,:,:,:)=Solver(x,y,x0-dx,y0,t,finalt,dx,dy,dt,dtK,Lpdf,m,n,K,r,M1,M2);
            flag(i-1,j)=1;
            cache(i,j,2)=MatTrajMetric(traj(i,j,:,:,:),traj(i-1,j,:,:,:),method,decoder);
        end
        if j>1 && flag(i,j-1)==0
            Updf = (1 / (4 * pi * D0 * dt0) * exp(-1 / (4 * D0 * dt0) * ((x - x0).^2 + (y - y0+dy).^2)));
            traj(i,j-1,:,:,:)=Solver(x,y,x0,y0-dy,t,finalt,dx,dy,dt,dtK,Updf,m,n,K,r,M1,M2);
            flag(i,j-1)=1;
            cache(i,j,1)=MatTrajMetric(traj(i,j,:,:,:),traj(i,j-1,:,:,:),method,decoder);
        end
        if i<m
            Rpdf = (1 / (4 * pi * D0 * dt0) * exp(-1 / (4 * D0 * dt0) * ((x - x0-dx).^2 + (y - y0).^2)));
            traj(i+1,j,:,:,:)=Solver(x,y,x0+dx,y0,t,finalt,dx,dy,dt,dtK,Rpdf,m,n,K,r,M1,M2);
            flag(i+1,j)=1;
            right=MatTrajMetric(traj(i,j,:,:,:),traj(i+1,j,:,:,:),method,decoder);
        end
        if j<n
            Dpdf = (1 / (4 * pi * D0 * dt0) * exp(-1 / (4 * D0 * dt0) * ((x - x0).^2 + (y - y0-dy).^2)));
            traj(i,j+1,:,:,:)=Solver(x,y,x0,y0+dy,t,finalt,dx,dy,dt,dtK,Dpdf,m,n,K,r,M1,M2);
            flag(i,j+1)=1;
            down=MatTrajMetric(traj(i,j,:,:,:),traj(i,j+1,:,:,:),method,decoder);
        end
        sol(i,j)=max(max(cache(i,j,1),down),max(cache(i,j,2),right));
        cache(i+1,j,2)=right;
        cache(i,j+1,1)=down;
    end
    if id_x*m_p+i>m
        break
    end
end

runtime = (cputime - starttime) / 60

%figure()
%contourf(x,y,log(sol)/finalt);

% [phiy, phix] = gradient(PHI, dy, dx);
% [psiy, psix] = gradient(PSI, dy, dx);

% %     phix=extendboundary(phix);
% %     psix=extendboundary(psix);
% %     phiy=extendboundary(phiy);
% %     psiy=extendboundary(psiy);

% a1 = phix.^2 + psix.^2;
% b1 = phix .* phiy + psix .* psiy;
% c1 = phiy.^2 + psiy.^2;
% lambda1 = (a1 + c1 + sqrt((a1 + c1).^2 - 4 * (a1 .* c1 - b1.^2))) / 2;
% FTLEd = log(max(1e-30, lambda1)) / 2 / abs(finalt);

clear traj;

[str]=strcat("./ParallelSol/DoubleGyre_MatTVNormParaCompuHighRank129x64D0001_",int2str(id),".mat");

save(str);
