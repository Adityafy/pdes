function [psi, omz, zeta, ...
                u, v, dpsi1, domz1, dzeta1, ...
                    fpv, fpvmag, lam1inst, lam1, res] ...
                        = gshFpvFOEAtIntervals(p,psitr,omztr,zetatr)
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
dt = p.dt_ts;
nmax = p.tstu/dt;
TsIntervals = p.TsIntervals;
N = Nx*Ny;

rng(p.seed,"twister");

%% Initial conditions and preallocation
% Running the transients
% [psitr, omztr, zetatr, ~, ~] = gshtric(p);

psimat = psitr(:,:,end);
zetamat = zetatr(:,:,end);
omzmat = omztr(:,:,end);

%% Initiating tangent space variables
% [dpsi1mat,domz1mat,fpv] = gshfpvic(p);


[dpsimat, ~, domzmat, ~, ~, ~, ~, ~, ~] = tsics(p);

dpsi1 = dpsimat(:,:,1);
domz1 = domzmat(:,:,1);

dpsi1mat = dpsimat(:,:,1);
domz1mat = domzmat(:,:,1);

fpv = [latToVec(dpsi1(:,:,1)); latToVec(domz1(:,:,1))];
fpvmag = norm(fpv);

dzeta1mat = zeros(Nx,Ny);

tN = p.tN;
nnorm = tN/dt;
lam1inst = [];
% res = zeros(1, tmax);

fprintf('\nCalculating coefficient matrices for implicit calcs...\n');
[matdivpsi, matdivomz] = impMatGSH(p,p.dt_ts);
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

    dpsi1mat = dpsi1mat + dt * ts1exp(psimat,zetamat, ...
                                    dpsi1mat,dzeta1mat,p,@sfdd);

    domz1mat = domz1mat + dt * ts2exp(psimat,domz1mat,dpsi1mat,p,@sfdd);

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
    
    if rem(n,nmax/TsIntervals) == 0
        %------------- saving dynamics at intervals ---------------
        psi(:,:,round(n*TsIntervals/nmax)) = psimat;
        omz(:,:,round(n*TsIntervals/nmax)) = omzmat;
        zeta(:,:,round(n*TsIntervals/nmax)) = zetamat;
        u(:,:,round(n*TsIntervals/nmax)) = ulat;
        v(:,:,round(n*TsIntervals/nmax)) = vlat;
        %----------- saving TS dynamics at intervals --------------
        dpsi1(:,:,round(n*TsIntervals/nmax)) = dpsi1mat;
        domz1(:,:,round(n*TsIntervals/nmax)) = domz1mat;
        dzeta1(:,:,round(n*TsIntervals/nmax)) = dzeta1mat;
    end

    if rem(n,nmax/10) == 0
        fprintf('This is time step %i, ',n);
        toc
    end
    % imagesc(dpsi1mat); colorbar; colormap jet; axis square;
    % plot(fpvmag(1:n+1),'-o');
    % drawnow;
end

lam1 = mean(lam1inst);

end