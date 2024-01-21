% clear all;
close all;
format compact
format long

m = 257;
D0 = 0.01;

n = (m + 1) / 2;
xmin = 0;
xmax = 2;
ymin = 0;
ymax = 1;
finalt = 5;

[x, y, dx, dy] = initial(m, n, xmin, xmax, ymin, ymax);
dt = 0.05;

miu = D0 * dt / dx^2;
M1 = left_mat(miu, m, n);
M2 = right_mat(miu, m, n);

i = (m + 3) / 4;
j = (3 * n + 1) / 4;
x0 = x(i, 1);
y0 = y(1, j);

velo=velocity(1);

%%%%%%%%\Phi_{t_0}^{t_0+0.5\Delta t}(x_0,y_0)%%%%%%%%%%%%
[u, v] = velo(x, y, 0, finalt);
[u1, v1] = velo(x, y, 0.5 * dt, finalt);
[xm1, xp1, ym1, yp1] = preparediff(u1);
[up, um, vp, vm] = WENO2(u1, xm1, xp1, ym1, yp1, dx, dy);
phi = x + 0.5 * 0.5 * dt * (u + u1) + 0.5 * (0.5 * dt)^2 * (u .* up + v .* vp);
[xm1, xp1, ym1, yp1] = preparediff(v1);
[up, um, vp, vm] = WENO2(v1, xm1, xp1, ym1, yp1, dx, dy);
psi = y + 0.5 * 0.5 * dt * (v + v1) + 0.5 * (0.5 * dt)^2 * (u .* up + v .* vp);
phi_point = phi(i, j);
psi_point = psi(i, j);
%%%%%%%%\Phi_{t_0}^{t_0+0.5\Delta t}(x_0,y_0)%%%%%%%%%%%%%%%%

%%%%%%%%%%%\Phi_{t_0+\Delta t}^{t_0+0.5\Delta t} (x,y)%%%%%%%%%%%%%%%%
[u, v] = velo(x, y, 0.5 * dt, finalt);
[u1, v1] = velo(x, y, dt, finalt);
[vp, up] = gradient(u, dy, dx);
phi = x - 0.5 * (0.5 * dt) * (u + u1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
[vp, up] = gradient(v, dy, dx);
psi = y - 0.5 * (0.5 * dt) * (v + v1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
%%%%%%%%%%%\Phi_{t_0+\Delta t}^{t_0+0.5\Delta t} (x,y)%%%%%%%%%%%%%%%%

pdf = 1 / (4 * pi * D0 * dt) * exp(-1 / (4 * D0 * dt) * ((phi - phi_point).^2 + (psi - psi_point).^2));
%       pdf=normalize(pdf,dx,dy);

t = dt;
k = 1;

while (t < finalt)

    [t finalt]

    %               if(t+dt>finalt)
    %                    dt=finalt-t;
    %                   miu=D0*dt/dx^2;
    %                  M1=left_mat(miu,m,n);
    %                   M2=right_mat(miu,m,n);
    %              end
    %
    %%%%%%%%%%% step 1 %%%%%%%%%%%%%%%%
    [u, v] = velo(x, y, t, finalt);
    [u1, v1] = velo(x, y, t + dt / 2, finalt);
    [vp, up] = gradient(u, dy, dx);
    phi = x - 0.5 * (0.5 * dt) * (u + u1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
    [vp, up] = gradient(v, dy, dx);
    psi = y - 0.5 * (0.5 * dt) * (v + v1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
    pdf = interp2(x', y', pdf', phi', psi', 'cubic')';
    %               pdf=normalize(pdf,dx,dy);
    %%%%%%%%%%% step 1 %%%%%%%%%%%%%%%%

    %%%%%%%%%%%%% step 2 %%%%%%%%%%%%%%%%%%%%
    pdf = pdf';
    Vright = M2 * pdf(:);
    pdf1 = M1 \ Vright;
    pdf1 = reshape(pdf1, n - 2, m - 2);
    pdf1 = pdf1';
    pdf1 = [zeros(1, n - 2); pdf1; zeros(1, n - 2)];
    pdf = [zeros(m, 1), pdf1, zeros(m, 1)];
    pdf = extendboundary(pdf);
    %              pdf=normalize(pdf,dx,dy);
    %%%%%%%%%%%%% step 2 %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%% step 3 %%%%%%%%%%%%%%%%
    u = u1; v = v1;
    [u1, v1] = velo(x, y, t + dt, finalt);
    [xm1, xp1, ym1, yp1] = preparediff(u);
    [up, um, vp, vm] = WENO2(u, xm1, xp1, ym1, yp1, dx, dy);
    phi = x - 0.5 * (0.5 * dt) * (u + u1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
    [xm1, xp1, ym1, yp1] = preparediff(v);
    [up, um, vp, vm] = WENO2(v, xm1, xp1, ym1, yp1, dx, dy);
    psi = y - 0.5 * (0.5 * dt) * (v + v1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
    pdf = interp2(x', y', pdf', phi', psi', 'cubic')';
    %               pdf=normalize(pdf,dx,dy);
    %%%%%%%%%%% step 3 %%%%%%%%%%%%%%%%

    %
    %             if (floor(t)~=floor(t+dt))
    %             imagesc(x(:,1),y(1,:),pdf',[0 max(max(pdf))])
    %             colormap('jet')
    %             axis equal
    %             axis xy
    %             axis([0 2 0 1])
    %             colorbar
    %             filename=strcat('PDF_',num2str(m),'_',num2str(n),'_',num2str(D0*10000),'_t',num2str(floor(t+dt)));
    %             print('-dpsc2',filename)
    %             print('-djpeg',filename)
    %             end
    k = k + 1;
    t = k * dt;

end

%          save m129.mat
%  xx=[1/32,1/64,1/128,1/256,1/512];
% L1=[0.008605761642389*3.88,0.008605761642389,0.002201537348495,5.578373233138480e-04,1.215163537308923e-04];
