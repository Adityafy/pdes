function [dpsinp1, res] = rk2tsgsh1(psi,zeta,dpsi,dzeta,p,tsdynfunc,sfddfunc)
% dpsinp1 = rk2tsgsh1(psi,zeta,dpsi,dzeta,p,tsdynfunc,fdsfunc)
% np1 means n+1.
    dt = p.dt_fpv;
    k1 = tsdynfunc(psi,zeta,dpsi,dzeta,p,sfddfunc);
    %u1 = u+k1*dt;
    k2 = tsdynfunc(psi,zeta,dpsi+k1*dt,dzeta,p,sfddfunc);
    dpsinp1 = dpsi + dt*((k1+k2)/2);

    dpsidt = dpsinp1 - dpsi;
    dpsidt = dpsidt/dt;
    rhs = (k1+k2)./2;
    res = max(max(dpsidt - rhs));
end