#! octave -qf
clear all
close all
format compact
format long
%pkg load ltfat

starttime = cputime;

arglist=argv();
id = str2num(argv{1})

autoload('SFsolver','./ImprovedSolver/SFsolver.oct');
addpath(genpath('./OSCAR files/'));

%parameter
m = 301;
n = 151;
D0 = 0.0001;
adap_depth = 1;

%domain info
xmin = 0;
xmax = 50;
ymin = 0;
ymax = 25;
finalt = 50;
K=51;

%parallel frame parameter
m_p = 76;
n_p = 76;
id_x=mod(id,4);
id_y=floor(id/4);

%rank approximation
r=25;
decoder = @(a,t) {LowRankDecoder(squeeze(a{1})(t,:,:),r,m,n),a{2:end}};

%initialization
[x, y, dx, dy] = initial(m, n, xmin, xmax, ymin, ymax);
dt0 = 0.5;
dt = 0.5;
dtK = (finalt-dt0)/(K-1);
miu = D0 * dt / dx^2;
M1 = left_mat(miu, m, n);
M2 = right_mat(miu, m, n);
traj = zeros(n_p+2,K,r,m+n+1);
cache = zeros(n_p+1,1); #1 - n_p: left-cache, n_p+1: up-cache
sol = zeros(m_p,n_p);

%metric
method = @(u) @(d1,d2) DisplacedSlicedWasserstein(d1,d2,x,y,u,15);
%method = @(u) @(d1,d2) DisplacedTotalVariation(d1,d2,u);

%folder
folder_str="./ParallelSol/ToPaper/OSCARD0001/";

for i = 0:m_p+1
    for j = 0:n_p+1
        if i+j == 0
            continue
        end
        [id, i, j]
        
        if id_x*m_p+i>0 && id_x*m_p+i<=m && id_y*n_p+j>0 && id_y*n_p+j<=n
            x0 = x(id_x*m_p+i, 1);
            y0 = y(1, id_y*n_p+j);
            t = dt0;
            pdf = 1 / (4 * pi * D0 * dt0) * exp(-1 / (4 * D0 * dt0) * ((x - x0).^2 + (y - y0).^2));
            curr = SFsolver(x,y,t,finalt,dx,dy,dt,pdf,K,r,M1,M2);
        else
            curr = zeros(K,r,m+n+1);
        end
        down=0;
        right=0;
        if j==0 && i>1 && (id_x<3 || (id_x==3&&id_x*m_p+i<=m+1))  %zero-th row, cache bottom value
            if id_y == 0
                cache(n_p+1) = 0;
            else
                cache(n_p+1) = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{traj(j+2,:,:,:),id_x*m_p+i-1,id_y*n_p+j+1},K,method(2),decoder);
            end
        elseif j>0 && j<n_p+1 %non-zero-th or n_p+1-th row
            if i==1 %first col, not on the boundary, lack of right info, only cache value
                if id_x == 0
                    cache(j) = 0;
                else
                    cache(j) = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{curr,id_x*m_p+i,id_y*n_p+j},K,method(1),decoder);
                end
            elseif  i>1 && (id_x<3 || (id_x==3&&id_x*m_p+i<=m+1))
                right = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{curr,id_x*m_p+i,id_y*n_p+j},K,method(1),decoder);
                down = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{traj(j+2,:,:,:),id_x*m_p+i-1,id_y*n_p+j+1},K,method(2),decoder);
                sol(i-1,j) = max(max(cache(n_p+1),down),max(cache(j),right));
                cache(j) = right;
                cache(n_p+1) = down;
            end
        end
        traj(j+1,:,:,:) = curr;
    end
end

runtime = (cputime - starttime) / 60;

%clear up memory
clear traj;
clear decoder;
clear method;

%cache
[str]=strcat(folder_str,"Cache_",int2str(id),".mat");
save(str);

[str]=strcat(folder_str,"finalsol_",int2str(id),".mat");
save("-mat7-binary",str,"runtime","sol");
