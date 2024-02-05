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
method = @(u) @(d1,d2) DisplacedTotalVariation(d1,d2,u);

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
r=25;
decoder = @(M) LowRankDecoder(M,r,m,n);

%initialization
[x, y, dx, dy] = initial(m, n, xmin, xmax, ymin, ymax);
dt0 = 0.05;
dt = 0.05;
dtK = (finalt-dt0)/(K-1);
miu = D0 * dt / dx^2;
M1 = left_mat(miu, m, n);
M2 = right_mat(miu, m, n);
traj = zeros(n_p+2,K,r,m+n+1);
cache = zeros(n_p+1,1); #1 - n_p: left-cache, n_p+1: up-cache
sol = zeros(m_p,n_p);

%main loop
for i = 0:m_p+1
    for j = 0:n_p+1
        [i, j]
        if id_x*m_p+i>0 && id_x*m_p+i<=m && id_y*n_p+j>0 && id_y*n_p+j<=n
            x0 = x(id_x*m_p+i, 1);
            y0 = y(1, id_y*n_p+j);
            t = dt0;
            pdf = 1 / (4 * pi * D0 * dt0) * exp(-1 / (4 * D0 * dt0) * ((x - x0).^2 + (y - y0).^2));
            curr = Solver(x,y,x0,y0,t,finalt,dx,dy,dt,dtK,pdf,m,n,K,r,M1,M2);
        else
            curr = zeros(K,r,m+n+1);
        end
        down=0;
        right=0;
        if j==0 && i>0
            cache(n_p+1) = MatTrajMetric(traj(j+1,:,:,:),traj(j+2,:,:,:),method(2),decoder);
        elseif j<n_p+1
            if i==1
                cache(j) = MatTrajMetric(traj(j+1,:,:,:),curr,method(1),decoder);
            elseif i>0
                right = MatTrajMetric(traj(j+1,:,:,:),curr,method(1),decoder);
                down = MatTrajMetric(traj(j+1,:,:,:),traj(j+2,:,:,:),method(2),decoder);
                sol(i-1,j) = max(max(cache(n_p+1),down),max(cache(j),right));
                cache(j) = right;
                cache(n_p+1) = down;
            end
        end
        traj(j+1,:,:,:) = curr;
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

[str]=strcat("./ParallelSol/DoubleGyre_MatTVnormParaCompu129x64D0001HighRank_",int2str(id),".mat");

save(str);
