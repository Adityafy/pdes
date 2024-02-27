function dpsirhs = ts1exp(psi,zeta,dpsi,dzeta,p,sfddfunc)
% dpsirhs = ts1exp(psi,zeta,dpsi,dzeta,p,fdsfunc)
% Calculates the RHS of the first tangent space equation for GSH
% calculations, via the explicit time marching methods. 
% Dependent on the spatial finite difference derivative.
    Nx = p.Nx;
    Ny = p.Ny;

    eps = p.eps;

    I = p.I;
    J = p.J;

    dpsirhs = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            dzetady = sfddfunc(zeta,0,1,i,j,p);
            dzetadx = sfddfunc(zeta,1,0,i,j,p);
            dpsidy = sfddfunc(psi,0,1,i,j,p);
            dpsidx = sfddfunc(psi,1,0,i,j,p);

            ddzetady = sfddfunc(dzeta,0,1,i,j,p);
            ddzetadx = sfddfunc(dzeta,1,0,i,j,p);
            ddpsidy = sfddfunc(dpsi,0,1,i,j,p);
            ddpsidx = sfddfunc(dpsi,1,0,i,j,p);
            
            d4dpsidx4 = sfddfunc(dpsi,4,0,i,j,p);
            d4dpsidx2dy2 = sfddfunc(dpsi,2,2,i,j,p);
            d4dpsidy4 = sfddfunc(dpsi,0,4,i,j,p);

            d2dpsidx2 = sfddfunc(dpsi,2,0,i,j,p);
            d2dpsidy2 = sfddfunc(dpsi,0,2,i,j,p);

            dpsirhs(i,j) = - dzetady * ddpsidx + dzetadx * ddpsidy ...
                       - ddzetady * dpsidx + ddzetadx * dpsidy ...
                       - d4dpsidx4 - 2 * d4dpsidx2dy2 - d4dpsidy4 ...
                       - 2 * d2dpsidx2 - 2 * d2dpsidy2 ...
                       + (eps - 1 - 3 * psi(I(i),J(j))^2) * dpsi(I(i),J(j));
        end
    end
end