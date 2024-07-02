function [pdf_o] = SubSolver(pdf,x,y,dx,dy,dt,M1,M2,u,v,u1,v1,u2,v2)
    [m,n] = size(pdf);
    [vp, up] = gradient(u, dy, dx);
    phi = x - 0.5 * (0.5 * dt) * (u + u1) +  0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
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
    u1 = u2;
    v1 = v2;
    [xm1, xp1, ym1, yp1] = preparediff(u);
    [up, um, vp, vm] = WENO2(u, xm1, xp1, ym1, yp1, dx, dy);
    %              phi=x-0.5*dt*(u+u1)+0.5*dt^2*(u1.*up+v1.*vp);
    phi = x - 0.5 * (0.5 * dt) * (u + u1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
    [xm1, xp1, ym1, yp1] = preparediff(v);
    [up, um, vp, vm] = WENO2(v, xm1, xp1, ym1, yp1, dx, dy);
    psi = y - 0.5 * (0.5 * dt) * (v + v1) + 0.5 * (0.5 * dt)^2 * (u1 .* up + v1 .* vp);
    pdf = interp2(x', y', pdf', phi', psi', 'cubic',0)';
    pdf = PDFnormalize(pdf, dx, dy);
    pdf_o = pdf;