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
    if id_x==0||id_x==1
        [str]=strcat("./ParallelSol/DoubleGyre_MatDisSqrtWeightTVnormParaCompu257*128D000001HighRank_",int2str(id),".mat");
        load(str);
        sol_com(id_x*m_d+1:(id_x+1)*m_d,id_y*n_d+1:(id_y+1)*n_d)=sol;
    else
        sol_com(id_x*m_d+1:(id_x+1)*m_d,id_y*n_d+1:(id_y+1)*n_d)=zeros(m_d,n_d);
    end
    %m_d=33;
end
sol_com=sol_com(5:m_d*2,5:n-5);
figure()
contourf(x(5:m_d*2,5:n-5),y(5:m_d*2,5:n-5),log(sol_com/dx)/10);
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