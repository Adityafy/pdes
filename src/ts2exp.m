function domzrhs = ts2exp(~,psi,domz,dpsi,p,sfddfunc)
% domzrhs = ts2exp(~,psi,domz,dpsi,p,fdsfunc)
% Calculates the RHS of the second tangent space equation for GSH
% calculations, via the explicit time marching methods. 
% Dependent on the spatial finite difference derivative.
    Nx = p.rmesh.Nx;
    Ny = p.rmesh.Ny;
    
    sig = p.con.sig;
    c = p.con.c;
    gm = p.con.gm;

    I = p.idx.I;
    J = p.idx.J;

    domzrhs = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            
            d2domzdx2 = sfddfunc(domz,2,0,i,j,p);
            d2domzdy2 = sfddfunc(domz,0,2,i,j,p);

            dpsidy = sfddfunc(psi,0,1,i,j,p);
            dpsidx = sfddfunc(psi,1,0,i,j,p);

            ddpsidy = sfddfunc(dpsi,0,1,i,j,p);
            ddpsidx = sfddfunc(dpsi,1,0,i,j,p);
            
            d3dpsidx3 = sfddfunc(dpsi,3,0,i,j,p);
            d3dpsidxdy2 = sfddfunc(dpsi,1,2,i,j,p);
            d3dpsidydx2 = sfddfunc(dpsi,2,1,i,j,p);
            d3dpsidy3 = sfddfunc(dpsi,0,3,i,j,p);
            
            d3psidx3 = sfddfunc(psi,3,0,i,j,p);
            d3psidxdy2 = sfddfunc(psi,1,2,i,j,p);
            d3psidydx2 = sfddfunc(psi,2,1,i,j,p);
            d3psidy3 = sfddfunc(psi,0,3,i,j,p);

            domzrhs(i,j) = sig * d2domzdx2 + sig * d2domzdy2 ...
                        - (c^2) * domz(I(i),J(j)) ...
                       - gm * dpsidy * d3dpsidx3 - gm * dpsidy * d3dpsidxdy2 ...
                       + gm * dpsidx * d3dpsidydx2 + gm * dpsidx * d3dpsidy3 ...
                       - gm * ddpsidy * d3psidx3 - gm * ddpsidy * d3psidxdy2 ...
                       + gm * ddpsidx * d3psidydx2 + gm * ddpsidx * d3psidy3 ;
        end
    end
end