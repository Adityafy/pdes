function domzrhs = ts2exp(~,psi,domz,dpsi,p,fdsfunc)
% domzrhs = ts2exp(~,psi,domz,dpsi,p,fdsfunc)
% Calculates the RHS of the first tangent space equation for GSH
% calculations, for the explicit time marching methods. 
% Dependent on the spatial finite difference derivative.
    Nx = p.Nx;
    Ny = p.Ny;
    
    sig = p.sig;
    c = p.c;
    gm = p.gm;

    I = p.I;
    J = p.J;

    domzrhs = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            
            d2domzdx2 = fdsfunc(domz,2,0,i,j,p);
            d2domzdy2 = fdsfunc(domz,0,2,i,j,p);

            dpsidy = fdsfunc(psi,0,1,i,j,p);
            dpsidx = fdsfunc(psi,1,0,i,j,p);

            ddpsidy = fdsfunc(dpsi,0,1,i,j,p);
            ddpsidx = fdsfunc(dpsi,1,0,i,j,p);
            
            d3dpsidx3 = fdsfunc(dpsi,3,0,i,j,p);
            d3dpsidxdy2 = fdsfunc(dpsi,1,2,i,j,p);
            d3dpsidydx2 = fdsfunc(dpsi,2,1,i,j,p);
            d3dpsidy3 = fdsfunc(dpsi,0,3,i,j,p);
            
            d3psidx3 = fdsfunc(psi,3,0,i,j,p);
            d3psidxdy2 = fdsfunc(psi,1,2,i,j,p);
            d3psidydx2 = fdsfunc(psi,2,1,i,j,p);
            d3psidy3 = fdsfunc(psi,0,3,i,j,p);

            domzrhs(i,j) = sig * d2domzdx2 + sig * d2domzdy2 ...
                        - (c^2) * domz(I(i),J(j)) ...
                       - gm * dpsidy * d3dpsidx3 - gm * dpsidy * d3dpsidxdy2 ...
                       + gm * dpsidx * d3dpsidydx2 + gm * dpsidx * d3dpsidy3 ...
                       - gm * ddpsidy * d3psidx3 - gm * ddpsidy * d3psidxdy2 ...
                       + gm * ddpsidx * d3psidydx2 + gm * ddpsidx * d3psidy3 ;
        end
    end
end