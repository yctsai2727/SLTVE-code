#! octave -qf
clear all
close all
format compact
format long

starttime = cputime;

arglist=argv();
id = str2num(argv{1})

%folder
folder_str_ad = "./ParallelSol/ToPaper/WassersteinD01hires/";

[str]=strcat(folder_str,"Cache_",int2str(id),".mat");
load(str);

sol_v=reshape(sol,m_p*n_p,1);
rho=prctile(sol_v,95);
Filter=double(sol>rho);

m_ad=2*m-1;
n_ad=2*n;
m_p_ad=2*m_p-1;
n_p_ad=2*n_p;
[x_ad, y_ad, dx,dy] = initial(m_ad, n_ad, xmin, xmax, ymin, ymax);

temp=zeros(m_p_ad,n_p_ad);
temp(1:2:m_p_ad,1:2:n_p_ad)=double(Filter);
for i=1:m_p_ad
    for j=1:n_p_ad
        if temp(i,j)==1
            if i>0
                temp(i-1,j)=2
            end
            if j>0
                temp(i,j-1)=2
            end
            if i<m_p_ad
                temp(i+1,j)=2
            end
            if j<n_p_ad
                temp(i,j+1)=2
            end
        end
    end
end
Filter = temp;

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
method = @(u) @(d1,d2) DisplacedSlicedWasserstein(d1,d2,x,y,u,15);
decoder = @(a,t) {LowRankDecoder(squeeze(a{1})(t,:,:),r,m,n),a{2:end}};

sol_ad = zeros(m_p,n_p);
sol_ad(1:2:m_p_ad,1:2:n_p_ad) = sol;
sol = sol_ad;

% if id_x==3
%     temp = interp2(x(id_x*m_p+1:m,id_y*n_p+1:(id_y+1)*n_p)',y(id_x*m_p+1:m,id_y*n_p+1:(id_y+1)*n_p)',sol(1:m_p-3,1:n_p)',x_ad(id_x*m_p_ad+1:m_ad,id_y*n_p_ad+1:(id_y+1)*n_p_ad)',y_ad(id_x*m_p_ad+1:m_ad,id_y*n_p_ad+1:(id_y+1)*n_p_ad)',"linear")';
%     sol=zeros(m_p_ad,n_p_ad);
%     sol(1:m_p_ad-3,1:n_p_ad) = temp;
% else
%     sol = interp2(x(id_x*m_p+1:(id_x+1)*m_p,id_y*n_p+1:(id_y+1)*n_p)',y(id_x*m_p+1:(id_x+1)*m_p,id_y*n_p+1:(id_y+1)*n_p)',sol',x_ad(id_x*m_p_ad+1:(id_x+1)*m_p_ad,id_y*n_p_ad+1:(id_y+1)*n_p_ad)',y_ad(id_x*m_p_ad+1:(id_x+1)*m_p_ad,id_y*n_p_ad+1:(id_y+1)*n_p_ad)',"linear")';
% end

for i = 1:m_p
    for j = 1:n_p
        [id, ad-1, i, j]
        
        if mod(i,2)==1 || mod(j,2)==1
            continue;
        end

        if Filter(i,j) != 1
            sol(i,j) = sum(sol)
        end

    end
end