m=129;
n=(m-1)/2;
m_d=33;
n_d=32;
sol_com=zeros(m_d*4,n_d*2);
[x, y, dx, dy] = initial(m, n, 0,2,0,1);
for id=0:7 
    id_x=mod(id,4);
    id_y=floor(id/4);
    % if id_x==3
    %     m_d=129-m_d*3;
    % end
    [str]=strcat("./ParallelSol/DoubleGyre_MatTVNormParaCompuHighRank129x64D0001_",int2str(id),".mat");
    load(str);
    % figure()
    % contourf(x,y,log(sol)/10);
    sol_com(id_x*m_d+1:(id_x+1)*m_d,id_y*n_d+1:(id_y+1)*n_d)=sol;
    %m_d=33;
end
sol_com=sol_com(10:m-10,10:n-10);
figure()
contourf(x(10:m-10,10:n-10),y(10:m-10,10:n-10),log(sol_com/dx)/10);
colorbar