function psitilde = rk2gsh1(psi,zeta,p,dynfunc)
% psitilde = rk2gsh1(psi,zeta,p,dynfunc)
    dt = p.dt;
    k1 = dynfunc(psi,zeta,p);
    %u1 = u+k1*dt;
    k2 = dynfunc(psi+k1*dt,zeta,p);
    psitilde = psi + dt*((k1+k2)/2);
end