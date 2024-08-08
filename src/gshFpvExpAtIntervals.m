function [psi, omz, zeta, ...
                u, v, dpsi1, domz1, dzeta1, ...
                    fpv, fpvmag, lam1inst, lam1, res] ...
                        = gshFpvExpAtIntervals(p,psitr,omztr,zetatr)
% [psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
%                   fpv, fpvmag, lam1inst, lam1, res] = gshFpvExpAtIntervals(p)
% Evolution of the first perturbation vector (fpv) that uses an explicit
% method. 

% run_type = 1; % this determines if the run is for transients (0) or FPV (1)

res_calc = 0; % put this = 1 if residual for rk2gsh1 has to be calculated
% therefore making res = 0 by default
res = 0;

% parameters
Nx = p.Nx;
Ny = p.Ny;
dt = p.dt_fpv;
nmax = p.totimeu/dt;
ToIntervals = p.ToIntervals;
N = Nx*Ny;

rng(p.seed,"twister");

%% Initial conditions and preallocation
% Running the transients
% [psitr, omztr, zetatr, ~, ~] = gshtric(p);

psimat = psitr(:,:,end);
zetamat = zetatr(:,:,end);
omzmat = omztr(:,:,end);

%% Initiating tangent space variables
[dpsi1mat,domz1mat,fpv] = gshfpvic(p);

fpvmag = norm(fpv);

dzeta1mat = zeros(Nx,Ny);

tN = p.tN;
nnorm = tN/dt;
lam1inst = [];
% res = zeros(1, tmax);

fprintf('\nCalculating coefficient matrices for implicit calcs...\n');
[matdivpsi, matdivomz] = impMatGSH(p,p.dt_fpv);
toc

fprintf('\nCalculating dynamics and the 1st Pert Vector...\n');

for n = 1:nmax
    %------------------ zeta, iterative ------------------------
    zetamat = iterativeZetaOmz(omzmat,zetamat,p);
    [ulat,vlat] = uvzeta(zetamat,p);
    
    %------------------ explicit nonlinear --------------------
    psitilde = rk2gsh1(psimat,zetamat,p,@nlpgsh1,dt); 
    psitilde = latToVec(psitilde);

    omztilde = rk2gsh2(omzmat,psimat,p,@nlpgsh2,dt);
    omztilde = latToVec(omztilde);
    
    %------------------ implicit CN linear ----------------------
    psivec = matdivpsi * psitilde;
    omzvec = matdivomz * omztilde;

    %--------- converting from vector to matrix ----------
    psimat = vecToLat(psivec,Nx,Ny);
    omzmat = vecToLat(omzvec,Nx,Ny);

    %------------------ explicit TS linear ----------------------
    dzeta1mat = iterativeZetaOmz(domz1mat,dzeta1mat,p);
    if res_calc == 0
        [dpsi1mat, ~] = rk2tsgsh1(psimat,zetamat,dpsi1mat, ...
                        dzeta1mat,p,@ts1exp,@sfdd);
    elseif res_calc == 1
        [dpsi1mat, res(n)] = rk2tsgsh1(psimat,zetamat,dpsi1mat, ...
                        dzeta1mat,p,@ts1exp,@sfdd);
    end
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
    
    if rem(n,nmax/ToIntervals) == 0
        %------------- saving dynamics at intervals ---------------
        psi(:,:,round(n*ToIntervals/nmax)) = psimat;
        omz(:,:,round(n*ToIntervals/nmax)) = omzmat;
        zeta(:,:,round(n*ToIntervals/nmax)) = zetamat;
        u(:,:,round(n*ToIntervals/nmax)) = ulat;
        v(:,:,round(n*ToIntervals/nmax)) = vlat;
        %----------- saving TS dynamics at intervals --------------
        dpsi1(:,:,round(n*ToIntervals/nmax)) = dpsi1mat;
        domz1(:,:,round(n*ToIntervals/nmax)) = domz1mat;
        dzeta1(:,:,round(n*ToIntervals/nmax)) = dzeta1mat;
    end

    if rem(n,nmax/10) == 0
        fprintf('This is time step %i, ',n);
        toc
    end
    
end

lam1 = mean(lam1inst);

end