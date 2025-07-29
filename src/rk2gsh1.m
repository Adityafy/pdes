function psitilde = rk2gsh1(psi,zeta,p,dynfunc,h)
% psitilde = rk2gsh1(psi,zeta,p,dynfunc,run_type)
k1 = dynfunc(psi,zeta,p);
k2 = dynfunc(psi+k1*h,zeta,p);
psitilde = psi + h*((k1+k2)/2);
end