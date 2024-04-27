#! octave -qf
clear all
close all
format compact
format long
#pkg load ltfat

starttime = cputime;

arglist=argv();
id = str2num(argv{1})

% autoload('c_emd','./pyemd/c_emd.mex');
% addpath(genpath('./pyemd/'));

%parameter
m = 129;
n = (m - 1) / 2;
D0 = 0.000001;
adap_depth = 1;

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
decoder = @(a,t) {LowRankDecoder(squeeze(a{1})(t,:,:),r,m,n),a{2:end}};

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
Filter=ones(m_p,n_p);

%metric
method = @(u) @(d1,d2) DisplacedSlicedWasserstein(d1,d2,x,y,u);

%folder
folder_str="./ParallelSol/SanityCheck/";

for ad=1:adap_depth
    %main loop
    % if ad==1
    %     [str]=strcat(folder_str,"Depth",int2str(ad-1),"Cache_",int2str(id),".mat");
    %     load(str);
    % else
    for i = 0:m_p+1
        for j = 0:n_p+1
            if i+j == 0
                continue
            end
            [id, ad-1, i, j]
            
            if id_x*m_p+i>0 && id_x*m_p+i<=m && id_y*n_p+j>0 && id_y*n_p+j<=n && (sum(Filter(max(1,i-1):min(i+1,m_p),min(max(j,1),n_p)))>=1 || sum(Filter(min(max(i,1),m_p),max(1,j-1):min(j+1,n_p)))>=1)
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
            if j==0 && i>1 && (id_x<3 || (id_x==3&&i<=m_p-2)) && Filter(i-1,j+1)==1  %zero-th row and point below is anchor, cache bottom value
                cache(n_p+1) = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{traj(j+2,:,:,:),id_x*m_p+i-1,id_y*n_p+j+1},K,method(2),decoder);
            elseif j>0 && j<n_p+1 %non-zero-th or n_p+1-th row
                if i==1 && id_x>0 && Filter(i,j)==1 %first col, not on the boundary, lack of right info, only cache value
                    cache(j) = MatTrajMetric({traj(j+1,:,:,:),id_x*m_p+i-1,id_y*n_p+j},{curr,id_x*m_p+i,id_y*n_p+j},K,method(1),decoder);
                elseif  i>1 && (id_x<3 || (id_x==3&&i<=m_p-2)) && Filter(i-1,j)==1
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
    
    %clear up memory
    clear traj;
    clear decoder;
    clear method;
    clear sample_v
    
    %cache
    [str]=strcat(folder_str,"Depth",int2str(ad-1),"Cache_",int2str(id),".mat");
    save(str);
    %end
    %Construct Filter
    sol_v=reshape(sol,m_p*n_p,1);
    rho=prctile(sol_v,95);
    Filter=double(sol>rho);

    %enlarge and interpolate current computation frame
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
    % for i=1:m_p     %remove negatively interpolated point, interpolate by only neighbouring point
    %     for j=1:n_p
    %         if (id_x<3||id_x==3&&i<=m_p-3)&&sol(i,j)<=0
    %             sol(i,j)=sol(max(i-1,1),j)+(id<3)*sol(min(i+1,m_p),j)+(id_x==3)*sol(min(i+1,m_p-3),j)+sol(i,max(j-1,1))+sol(i,min(j+1,n_p));
    %         end
    %     end
    % end

    %adjust filter to fit the new frame
    temp=zeros(m_p_ad,n_p_ad);
    temp(1:2:m_p_ad,1:2:n_p_ad)=Filter;
    for i=1:m_p_ad
        for j=1:n_p_ad
            if temp(i,j)==1
                if i>1
                    temp(i-1,j)=2;
                end
                if j>1
                    temp(i,j-1)=2;
                end
                if ((id_x<3&&i<m_p_ad) || (id_x==3&&i<m_p_ad-3)) && temp(i+1,j)<1
                    temp(i+1,j)=2;
                end
                
                if j<n_p && temp(i,j+1)<1
                    temp(i,j+1)=2;
                end
            end
        end
    end
    Filter=double(temp>=1);

    %update parameter
    x=x_ad;
    y=y_ad;
    m=m_ad;
    n=n_ad;
    m_p=m_p_ad;
    n_p=n_p_ad;
    miu = D0 * dt / dx^2;
    M1 = left_mat(miu, m, n);
    M2 = right_mat(miu, m, n);
    traj = zeros(n_p+2,K,r,m+n+1);
    cache = zeros(n_p+1,1);
    method = @(u) @(d1,d2)  femd(d1,d2);
    decoder = @(a,t) LowRankDecoder(squeeze(a{1})(t,:,:),r,m,n);
end

runtime = (cputime - starttime) / 60

clear traj;

[str]=strcat(folder_str,"finalsol_",int2str(id),".mat");

save("-mat7-binary",str,"runtime","sol");
