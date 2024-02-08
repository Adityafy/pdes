function omztilde = rk2gsh2(omz,psi,p,dynfunc)
% omztilde = rk2gsh2(omz,psi,p,dynfunc)
    dt = p.dt;
    k1 = dynfunc(omz,psi,p);
    %u1 = u+k1*dt;
    k2 = dynfunc(omz+k1*dt,psi,p);
    omztilde = omz + dt*((k1+k2)/2);
end