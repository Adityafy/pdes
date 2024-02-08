function domznp1 = rk2tsgsh2(omz,psi,domz,dpsi,p,tsdynfunc,fdsfunc)
% domznp1 = rk2tsgsh2(omz,psi,domz,dpsi,p,tsdynfunc,fdsfunc)
% np1 means n+1
    dt = p.dt;
    k1 = tsdynfunc(omz,psi,domz,dpsi,p,fdsfunc);
    %u1 = u+k1*dt;
    k2 = tsdynfunc(omz,psi,domz+k1*dt,dpsi,p,fdsfunc);
    domznp1 = domz + dt*((k1+k2)/2);
end