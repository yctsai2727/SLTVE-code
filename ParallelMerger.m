m=257;
n=(m-1)/2;
m_d=65;
n_d=64;
epsi=1e-10;
sol_com=zeros(m_d*4,n_d*2);
[x, y, dx, dy] = initial(m, n, 0,2,0,1);
for id=0:7 
    id_x=mod(id,4);
    id_y=floor(id/4);
    % if id_x==3
    %     m_d=129-m_d*3;
    % end
    [str]=strcat("./ParallelSol/DoubleGyre_MatDisTVnorm257x128adapD000001HighRank_",int2str(id),".mat");
    load(str);
    sol_com(id_x*m_d+1:(id_x+1)*m_d,id_y*n_d+1:(id_y+1)*n_d)=sol;
    %m_d=33;
end
sol_com=sol_com(10:m-10,10:n-10);
sol_com=sol_com+(sol_com<epsi)*epsi;
%sol_com.*imag(log(sol_com))>0
figure()
%contourf(x(10:m-10,10:n-10),y(10:m-10,10:n-10),log(sol_com/dx)/10);
imagesc(x(10:m-10,1),y(1,10:n-10),(log(sol_com/dx)/10)',[0,2]);
%contourf(x,y,sol_com);
colorbar

% id=3;
% id_x=mod(id,4);
% id_y=floor(id/4);
% [str]=strcat("./ParallelSol/DoubleGyre_MatDisWeightTVnormParaCompu129x64D0001HighRank_",int2str(3),".mat");
% load(str);
% figure()
% contourf(x(id_x*m_d+1:(id_x+1)*m_d-10,id_y*n_d+1:(id_y+1)*n_d),y(id_x*m_d+1:(id_x+1)*m_d-10,id_y*n_d+1:(id_y+1)*n_d),log(sol(1:m_d-10,1:n_d)/dx)/10);
% colorbar