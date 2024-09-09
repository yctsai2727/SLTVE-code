function [adap_sol] = Solver(sol,x,y,t0,finalt,dx,dy,dt,K,r)
[m,n] = size(x);
[m_p,n_p] = size(sol);
sol_v=reshape(sol,m_p*n_p,1);
rho=prctile(sol_v,95);
Filter=uint8(sol>rho);
