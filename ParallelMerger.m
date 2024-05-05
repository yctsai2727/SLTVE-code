m=65;
n=(m-1)/2;
m_d=17;
n_d=16;
epsi=1e-10;
sol_com=zeros(m_d*4,n_d*2);
[x, y, dx, dy] = initial(m, n, 0,2,0,1);
for id=0:7 
    id_x=mod(id,4);
    id_y=floor(id/4);
    % if id_x==3
    %     m_d=129-m_d*3;
    % end
    [str]=strcat("./ParallelSol/SlicedWasserstein/Depth0Cache_",int2str(id),".mat");
    load(str);
    sol_com(id_x*m_d+1:(id_x+1)*m_d,id_y*n_d+1:(id_y+1)*n_d)=sol;
    %m_d=33;
end
sol_com=sol_com(1:m,1:n);
%sol_com=sol_com.*(sol_com>epsi)+epsi*(sol_com<=epsi);
%sol_com.*imag(log(sol_com))>0
figure()
contourf(x(1:m,1:n),y(1:m,1:n),log(sol_com/dx)/10);
%imagesc(x(5:m-5,1),y(1,5:n-5),(log(sol_com/dx)/10)',[0,2]);
colorbar

% id=3;
% id_x=mod(id,4);
% id_y=floor(id/4);
% [str]=strcat("./ParallelSol/DoubleGyre_MatDisWeightTVnormParaCompu129x64D0001HighRank_",int2str(3),".mat");
% load(str);
% figure()
% contourf(x(id_x*m_d+1:(id_x+1)*m_d-10,id_y*n_d+1:(id_y+1)*n_d),y(id_x*m_d+1:(id_x+1)*m_d-10,id_y*n_d+1:(id_y+1)*n_d),log(sol(1:m_d-10,1:n_d)/dx)/10);
% colorbar