m=151;
n=76;
m_d=38;
n_d=38;
epsi=1e-10;
sol_com=zeros(m_d*4,n_d*2);
max_runtime = -1;
[x, y, dx, dy] = initial(m, n, 0,50,0,25);
for id=0:7 
    id_x=mod(id,4);
    id_y=floor(id/4);
    % if id_x==3
    %     m_d=129-m_d*3;
    % end
    if max_runtime < 0 || runtime > max_runtime
        max_runtime = runtime;
    end
    [str]=strcat("./ParallelSol/ToPaper/OSCARD0001/Cache_",int2str(id),".mat");
    load(str);
    sol_com(id_x*m_d+1:(id_x+1)*m_d,id_y*n_d+1:(id_y+1)*n_d)=sol;
    %m_d=33;
end
[max_runtime]
sol_com=sol_com(2:m-1,2:n-1);
%sum(sum(abs(imag(sol_com))>0))
%sol_com=sol_com.*(sol_com>epsi)+epsi*(sol_com<=epsi);
%sol_com.*imag(log(sol_com))>0
figure()
%contourf(x(2:m-1,2:n-1),y(2:m-1,2:n-1),log(sol_com/dx)/50);
imagesc(x(2:m-1,1),y(1,2:n-1),(log(sol_com/dx)/50)',[0,0.25]);
colorbar

% id=3;
% id_x=mod(id,4);
% id_y=floor(id/4);
% [str]=strcat("./ParallelSol/DoubleGyre_MatDisWeightTVnormParaCompu129x64D0001HighRank_",int2str(3),".mat");
% load(str);
% figure()
% contourf(x(id_x*m_d+1:(id_x+1)*m_d-10,id_y*n_d+1:(id_y+1)*n_d),y(id_x*m_d+1:(id_x+1)*m_d-10,id_y*n_d+1:(id_y+1)*n_d),log(sol(1:m_d-10,1:n_d)/dx)/10);
% colorbar