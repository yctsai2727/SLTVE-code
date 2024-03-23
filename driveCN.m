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
D0 = 0.000001;

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
r=16;
decoder = @(a,t) LowRankDecoder(squeeze(a{1})(t,:,:),r,m,n);

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

%metric
method = @(u) @(d1,d2) DisplacedTotalVariation(d1,d2,u);

%first main loop
% for i = 0:m_p+1
%     for j = 0:n_p+1
%         [id, i, j]
%         if id_x*m_p+i>0 && id_x*m_p+i<=m && id_y*n_p+j>0 && id_y*n_p+j<=n
%             x0 = x(id_x*m_p+i, 1);
%             y0 = y(1, id_y*n_p+j);
%             t = dt0;
%             pdf = 1 / (4 * pi * D0 * dt0) * exp(-1 / (4 * D0 * dt0) * ((x - x0).^2 + (y - y0).^2));
%             curr = Solver(x,y,x0,y0,t,finalt,dx,dy,dt,dtK,pdf,m,n,K,r,M1,M2);
%         else
%             curr = zeros(K,r,m+n+1);
%         end
%         down=0;
%         right=0;
%         if j==0 && i>1 && (id_x<3 || id_x==3&&i<=m_p-2) %first row, cache bottom value
%             cache(n_p+1) = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{traj(j+2,:,:,:),id_x*m_p+i-1,id_y*n_p+j+1},K,method(2),decoder);
%         elseif j>0 && j<n_p+1 %non-zero-th or n_p+1-th row
%             if i==1 %first col, lack of right info, only cache value
%                 cache(j) = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{curr,id_x*m_p+i,id_y*n_p+j},K,method(1),decoder);
%             elseif i>1 && (id_x<3 || id_x==3&&i<=m_p-2) %second+ col, compute LTVE of prev col
%                 right = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{curr,id_x*m_p+i,id_y*n_p+j},K,method(1),decoder);
%                 down = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{traj(j+2,:,:,:),id_x*m_p+i-1,id_y*n_p+j+1},K,method(2),decoder);
%                 sol(i-1,j) = max(max(cache(n_p+1),down),max(cache(j),right));
%                 cache(j) = right;
%                 cache(n_p+1) = down;
%             end
%         end
%         traj(j+1,:,:,:) = curr;
%     end
% end

clear traj;
clear decoder;
clear method;
[str]=strcat("./ParallelSol/TVNormTemp",int2str(id),".mat");
%save(str);

load(str);

%adaptive step
sol_v=reshape(sol(1:m_p-4,1:n_p),(m_p-4)*n_p,1);
rho=prctile(sol_v,85);
Filter=double(sol>rho);
for i=1:m_p
    for j=1:n_p
        if Filter(i,j)==1
            if i>1
                Filter(i-1,j)=2;
            end
            if (id_x<3&&i<m_p) || (id_x==3&&i<m_p-3)
                Filter(i+1,j)=2;
            end
            if j>1
                Filter(i,j-1)=2;
            end
            if j<n_p
                Filter(i,j+1)=2;
            end
        end
    end
end
Filter=Filter>=1;
m_ad=2*m-1;
n_ad=2*n;
m_p_ad=2*m_p-1;
n_p_ad=2*n_p;
[x_ad, y_ad, dx,dy] = initial(m_ad, n_ad, xmin, xmax, ymin, ymax);
if id_x==3
    temp = interp2(x(id_x*m_p+1:m,id_y*n_p+1:(id_y+1)*n_p)',y(id_x*m_p+1:m,id_y*n_p+1:(id_y+1)*n_p)',sol(1:m_p-3,1:n_p)',x_ad(id_x*m_p_ad+1:m_ad,id_y*n_p_ad+1:(id_y+1)*n_p_ad)',y_ad(id_x*m_p_ad+1:m_ad,id_y*n_p_ad+1:(id_y+1)*n_p_ad)',"spline")';
    sol=zeros(m_p_ad,n_p_ad);
    sol(1:m_p_ad-3,1:n_p_ad) = temp;
else
    sol = interp2(x(id_x*m_p+1:(id_x+1)*m_p,id_y*n_p+1:(id_y+1)*n_p)',y(id_x*m_p+1:(id_x+1)*m_p,id_y*n_p+1:(id_y+1)*n_p)',sol',x_ad(id_x*m_p_ad+1:(id_x+1)*m_p_ad,id_y*n_p_ad+1:(id_y+1)*n_p_ad)',y_ad(id_x*m_p_ad+1:(id_x+1)*m_p_ad,id_y*n_p_ad+1:(id_y+1)*n_p_ad)',"spline")';
end
%sol = sol.*(sol>=0);    #remove negatively interpolated point
for i=1:m_p     #remove negatively interpolated point, interpolate by only neighbor point
    for j=1:n_p
        if (id_x<3||id_x==3&&i<=m_p-3)&&sol(i,j)<=0
            sol(i,j)=sol(max(i-1,1),j)+(id<3)*(min(i+1,m_p),j)+(id_x==3)*(min(i+1,m_p-3),j)+sol(i,max(j-1,1))+sol(i,min(j+1,n_p));
        end
    end
end
x=x_ad;
y=y_ad;
m=m_ad;
n=n_ad;
m_p=m_p_ad;
n_p=n_p_ad;
temp=zeros(m_p,n_p);
temp(1:2:m_p,1:2:n_p)=Filter;
Filter=temp;
M1 = left_mat(miu, m, n);
M2 = right_mat(miu, m, n);
traj = zeros(n_p+2,K,r,m+n+1);
cache = zeros(n_p+1,1); #1 - n_p: left-cache, n_p+1: up-cache
method = @(u) @(d1,d2) DisplacedTotalVariation(d1,d2,u);
decoder = @(a,t) LowRankDecoder(squeeze(a{1})(t,:,:),r,m,n);

%second main loop
for i = 0:m_p+1
    for j = 0:n_p+1
        if i+j == 0
            continue
        end
        [i, j]
        
        if id_x*m_p+i>0 && id_x*m_p+i<=m && id_y*n_p+j>0 && id_y*n_p+j<=n && sum(Filter(max(1,i-1):min(i+1,m_p),max(1,j-1):min(j+1,n_p)))>=1
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
        if j==0 && i>1 && (id_x<3 || id_x==3&&i<=m_p-2) && Filter(i-1,j+1)==1  %zero-th row and point below is anchor, cache bottom value
            cache(n_p+1) = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{traj(j+2,:,:,:),id_x*m_p+i-1,id_y*n_p+j+1},K,method(2),decoder);
        elseif j>0 && j<n_p+1 %non-zero-th or n_p+1-th row
            if i==1 && Filter(i,j)==1 %first col, lack of right info, only cache value
                cache(j) = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{curr,id_x*m_p+i,id_y*n_p+j},K,method(1),decoder);
            elseif  i>1 && (id_x<3 || id_x==3&&i<=m_p-2) && Filter(i-1,j)==1
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


runtime = (cputime - starttime) / 60

clear traj;

[str]=strcat("./ParallelSol/DoubleGyre_MatDisTVnorm257x128adapD000001HighRank_",int2str(id),".mat");

save("-mat7-binary",str,"runtime","sol");
