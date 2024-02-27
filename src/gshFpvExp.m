function [dpsi1, domz1, dzeta1, fpv, fpvmag, lam1inst, lam1, res] = ...
gshFpvExp(p,psi,omz,zeta,evofunc1,evofunc2,fpvicfunc)
% [dpsi1, domz1, dzeta1, fpv, fpvmag, lam1inst, lam1] = gshFpvExp(p,psi,omz,zeta)
% Evolution of the first perturbation vector (fpv) that uses an explicit
% method. 
Nx = p.Nx;
Ny = p.Ny;
dt = p.dt;
tmax = p.tmax;
N = Nx*Ny;

rng(p.seed,"twister");

[dpsi1,domz1,fpv] = fpvicfunc(p);

fpvmag = norm(fpv(:,1));

dzeta1 = zeros(Nx,Ny,tmax);

normtime = p.normtime;
tnorm = normtime*dt;
lam1inst = [];
res = zeros(1, tmax);

fprintf('\nEvolving the first perturbation vector...\n');

for t = 1:tmax
    
    %------------------ explicit TS linear ----------------------
    dzeta1(:,:,t) = iterativeZetaOmz(domz1(:,:,t),dzeta1(:,:,t),p);
    [dpsi1(:,:,t+1), res(t)] = evofunc1(psi(:,:,t),zeta(:,:,t),dpsi1(:,:,t), ...
                        dzeta1(:,:,t),p,@ts1exp,@sfdd);
    domz1(:,:,t+1) = evofunc2(omz(:,:,t),psi(:,:,t),domz1(:,:,t), ...
                        dpsi1(:,:,t),p,@ts2exp,@sfdd);
    fpv(:,t+1) = [latToVec(dpsi1(:,:,t+1)); latToVec(domz1(:,:,t+1))];
    fpvmag(t+1) = norm(fpv(:,t+1));
    
    %------------------ renormalization and LLE ----------------------
    if rem(t,normtime) == 0
        lam1inst = [lam1inst, ...
            (1/tnorm) * log(abs(norm(fpv(:,t+1))))];
        fpv(:,t+1) = fpv(:,t+1)./norm(fpv(:,t+1));
        dpsi1(:,:,t+1) = vecToLat(fpv(1:N,t+1),Nx,Ny);
        domz1(:,:,t+1) = vecToLat(fpv(N+1:2*N,t+1),Nx,Ny);
    end
    
    dpsi1(:,:,t+1) = vecToLat(fpv(1:N,t+1),Nx,Ny);
    domz1(:,:,t+1) = vecToLat(fpv(N+1:2*N,t+1),Nx,Ny);

    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
    
end

lam1 = mean(lam1inst);

end