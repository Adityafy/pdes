function [psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
                    fpv, fpvmag, lam1inst, lam1, res] = gshFpvExpAtIntervals(p)
% [psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
%                   fpv, fpvmag, lam1inst, lam1, res] = gshFpvExpAtIntervals(p)
% Evolution of the first perturbation vector (fpv) that uses an explicit
% method. 
Nx = p.Nx;
Ny = p.Ny;
dt = p.dt;
tmax = p.tmax;
intervals = p.intervals;
N = Nx*Ny;

rng(p.seed,"twister");

%% Initial conditions and preallocation
% Running the transients
[psimatic, omzmatic, zetamatic, ~, ~] = gshtric(p);

psimat = psimatic;
zetamat = zetamatic;
omzmat = omzmatic;

%% Initiating tangent space variables
[dpsi1mat,domz1mat,fpv] = gshfpvic(p);

fpvmag = norm(fpv);

dzeta1mat = zeros(Nx,Ny);

tN = p.tN;
nnorm = tN/dt;
lam1inst = [];
% res = zeros(1, tmax);

fprintf('\nCalculating coefficient matrices for implicit calcs...\n');
[matdivpsi, matdivomz] = impMatGSH(p);
toc

fprintf('\nCalculating dynamics and the 1st Pert Vector...\n');

for n = 1:tmax
    %------------------ zeta, iterative ------------------------
    zetamat = iterativeZetaOmz(omzmat,zetamat,p);
    [ulat,vlat] = uvzeta(zetamat,p);
    
    %------------------ explicit nonlinear --------------------
    psitilde = rk2gsh1(psimat,zetamat,p,@nlpgsh1); 
    psitilde = latToVec(psitilde);

    omztilde = rk2gsh2(omzmat,psimat,p,@nlpgsh2);
    omztilde = latToVec(omztilde);
    
    %------------------ implicit CN linear ----------------------
    psivec = matdivpsi * psitilde;
    omzvec = matdivomz * omztilde;

    %--------- converting from vector to matrix ----------
    psimat = vecToLat(psivec,Nx,Ny);
    omzmat = vecToLat(omzvec,Nx,Ny);

    %------------------ explicit TS linear ----------------------
    dzeta1mat = iterativeZetaOmz(domz1mat,dzeta1mat,p);
    [dpsi1mat, res(n)] = rk2tsgsh1(psimat,zetamat,dpsi1mat, ...
                        dzeta1mat,p,@ts1exp,@sfdd);
    domz1mat = rk2tsgsh2(omzmat,psimat,domz1mat, ...
                        dpsi1mat,p,@ts2exp,@sfdd);
    fpv = [latToVec(dpsi1mat); latToVec(domz1mat)];
    fpvmag(n+1) = norm(fpv);
    
    %------------------ renormalization and LLE ----------------------
    if rem(n,nnorm) == 0
        lam1inst = [lam1inst, (1/tN) * log(abs(norm(fpv)))];
        fpv = fpv./norm(fpv);
        dpsi1mat = vecToLat(fpv(1:N),Nx,Ny);
        domz1mat = vecToLat(fpv(N+1:2*N),Nx,Ny);
    end
    
    dpsi1mat = vecToLat(fpv(1:N),Nx,Ny);
    domz1mat = vecToLat(fpv(N+1:2*N),Nx,Ny);
    
    if rem(n,tmax/intervals) == 0
        %------------- saving dynamics at intervals ---------------
        psi(:,:,round(n*intervals/tmax)) = psimat;
        omz(:,:,round(n*intervals/tmax)) = omzmat;
        zeta(:,:,round(n*intervals/tmax)) = zetamat;
        u(:,:,round(n*intervals/tmax)) = ulat;
        v(:,:,round(n*intervals/tmax)) = vlat;
        %----------- saving TS dynamics at intervals --------------
        dpsi1(:,:,round(n*intervals/tmax)) = dpsi1mat;
        domz1(:,:,round(n*intervals/tmax)) = domz1mat;
        dzeta1(:,:,round(n*intervals/tmax)) = dzeta1mat;
    end

    if rem(n,tmax/10) == 0
        fprintf('This is time step %i, ',n);
        toc
    end
    
end

lam1 = mean(lam1inst);

end