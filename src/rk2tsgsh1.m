function dpsinp1 = rk2tsgsh1(psi,zeta,dpsi,dzeta,p,tsdynfunc,fdsfunc)
% dpsinp1 = rk2tsgsh1(psi,zeta,dpsi,dzeta,p,tsdynfunc,fdsfunc)
% np1 means n+1.
    dt = p.dt;
    k1 = tsdynfunc(psi,zeta,dpsi,dzeta,p,fdsfunc);
    %u1 = u+k1*dt;
    k2 = tsdynfunc(psi,zeta,dpsi+k1*dt,dzeta,p,fdsfunc);
    dpsinp1 = dpsi + dt*((k1+k2)/2);
end