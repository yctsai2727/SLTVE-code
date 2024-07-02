function [sol] = Solver(x,y,t0,finalt,dx,dy,dt,pdf0,K,r,M1,M2)
    sub_start = cputime;
    [m,n] = size(x);
    flag = false;
    pdf=pdf0;
    k=2;
    c=2;
    sol=zeros(K,r,m+n+1);
    sol(1,:,:)=LowRankApprox(pdf0,r);
    jump = floor(floor((finalt-t0)/dt)/(K-2));
    t=t0;

    velo=OSCARwrapper(dt);

    while (t < finalt)
        %%%%%%%%%%%% step 1 %%%%%%%%%%%%%%%%
        [u, v] = velo(x, y, t, finalt);
        [u1, v1] = velo(x, y, t + dt / 2, finalt);
        [vp, up] = gradient(u, dy, dx);
        phi = x - 0.5 * (0.5 * dt) * (u + u1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
        [vp, up] = gradient(v, dy, dx);
        psi = y - 0.5 * (0.5 * dt) * (v + v1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
        pdf = interp2(x', y', pdf', phi', psi', 'cubic',0)';
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
        pdf = interp2(x', y', pdf', phi', psi', 'cubic',0)';
        pdf = PDFnormalize(pdf, dx, dy);
        %%%%%%%%%%% step 3 %%%%%%%%%%%%%%%%
        t = t + dt;
        if t>=finalt && !flag
            t = t-dt;
            dt = finalt - t;
            t = finalt;
            flag = true;
        end

        if mod(k-1,jump)==0 && c<K
            %sol=sol+sum(sum(abs(pdf))); %TV-norm
            sol(c,:,:)=LowRankApprox(pdf,r);
            %t-dt
            % if c==25 || c== 50 || c==75
            %     max(max(abs(pdf-LowRankDecoder(sol(c,:,:),r,m,n))))
            % end
            c=c+1;
            % sol=sol+Ex*Ex+Ey*Ey;
        end
        k=k+1;
    end
    sol(K,:,:)=LowRankApprox(pdf,r);
    solver_duration = cputime() - sub_start
    %N=(k-(mod(k,2)==1))/2;
    %sol=sqrt(sol/N);