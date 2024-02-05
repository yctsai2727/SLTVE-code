function [sol] = Solver(x,y,x0,y0,t0,finalt,dx,dy,dt,dtK,pdf0,m,n,K,r,M1,M2)

    pdf=pdf0;
    k=1;
    c=2;
    sol=zeros(K,r,m+n+1);
    sol(1,:,:)=LowRankApprox(pdf0,r);
    jump=ceil(dtK/dt);
    t=t0;

    velo=velocity(1);

    while (t < finalt)
        %%%%%%%%%%%% step 1 %%%%%%%%%%%%%%%%
        [u, v] = velo(x, y, t, finalt);
        [u1, v1] = velo(x, y, t + dt / 2, finalt);
        [vp, up] = gradient(u, dy, dx);
        phi = x - 0.5 * (0.5 * dt) * (u + u1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
        [vp, up] = gradient(v, dy, dx);
        psi = y - 0.5 * (0.5 * dt) * (v + v1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
        pdf = interp2(x', y', pdf', phi', psi', 'cubic')';
        pdf = PDFnormalize(pdf, dx, dy);
        %%%%%%%%%%% step 1 %%%%%%%%%%%%%%%%

        %%%%%%%%%%%%% step 2 %%%%%%%%%%%%%%%%%%%%
        pdf = pdf';
        Vright = M2 * pdf(:);
        pdf1 = M1 \ Vright;
        pdf1 = reshape(pdf1, n - 2, m - 2);
        pdf1 = pdf1';
        pdf1 = [zeros(1, n - 2); pdf1; zeros(1, n - 2)];
        pdf = [zeros(m, 1), pdf1, zeros(m, 1)];
        %             pdf=normalize(pdf,dx,dy);
        pdf = extendboundary(pdf);
        %%%%%%%%%%%%% step 2 %%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%% step 3 %%%%%%%%%%%%%%%%
        u = u1; v = v1;
        [u1, v1] = velo(x, y, t + dt, finalt);
        [xm1, xp1, ym1, yp1] = preparediff(u);
        [up, um, vp, vm] = WENO2(u, xm1, xp1, ym1, yp1, dx, dy);
        %              phi=x-0.5*dt*(u+u1)+0.5*dt^2*(u1.*up+v1.*vp);
        phi = x - 0.5 * (0.5 * dt) * (u + u1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
        [xm1, xp1, ym1, yp1] = preparediff(v);
        [up, um, vp, vm] = WENO2(v, xm1, xp1, ym1, yp1, dx, dy);
        psi = y - 0.5 * (0.5 * dt) * (v + v1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
        pdf = interp2(x', y', pdf', phi', psi', 'cubic')';
        pdf = PDFnormalize(pdf, dx, dy);
        %%%%%%%%%%% step 3 %%%%%%%%%%%%%%%%
        t = t + dt;
        if mod(k,jump)==0
            %sol=sol+sum(sum(abs(pdf))); %TV-norm
            sol(c,:,:)=LowRankApprox(pdf,r);
            % if c==25 || c== 50 || c==75
            %     max(max(abs(pdf-LowRankDecoder(sol(c,:,:),r,m,n))))
            % end
            c=c+1;
            % sol=sol+Ex*Ex+Ey*Ey;
        end
        k=k+1;
    end
    sol(K,:,:)=LowRankApprox(pdf,r);
    %N=(k-(mod(k,2)==1))/2;
    %sol=sqrt(sol/N);