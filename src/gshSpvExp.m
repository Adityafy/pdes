function [dpsi2, domz2, dzeta2, spv, spvmag, lam2inst, lam2, res] = ...
gshSpvExp(p,psi,omz,zeta,fpv,evofunc1,evofunc2,fpvicfunc)
% [dpsi2, domz2, dzeta2, spv, spvmag, lam2inst, lam2, res] = ...
%        gshSpvExp(p,psi,omz,zeta,fpv,evofunc1,evofunc2,fpvicfunc)
% Evolution of the second perturbation vector (spv) that uses an explicit
% method. 
Nx = p.Nx;
Ny = p.Ny;
dt = p.dt;
tmax = p.tmax;
N = Nx*Ny;

rng(p.seed+1,"twister");

[dpsi2,domz2,spv] = fpvicfunc(p);

spvmag = norm(spv(:,1));

dzeta2 = zeros(Nx,Ny,tmax);

normtime = p.normtime;
tnorm = normtime*dt;
lam2inst = [];
res = zeros(1, tmax);

fprintf('\nEvolving the second perturbation vector...\n');

for t = 1:tmax
    
    %------------------ explicit TS linear ----------------------
    dzeta2(:,:,t) = iterativeZetaOmz(domz2(:,:,t),dzeta2(:,:,t),p);
    [dpsi2(:,:,t+1), res(t)] = evofunc1(psi(:,:,t),zeta(:,:,t),dpsi2(:,:,t), ...
                        dzeta2(:,:,t),p,@ts1exp,@sfdd);
    domz2(:,:,t+1) = evofunc2(omz(:,:,t),psi(:,:,t),domz2(:,:,t), ...
                        dpsi2(:,:,t),p,@ts2exp,@sfdd);
    spv(:,t+1) = [latToVec(dpsi2(:,:,t+1)); latToVec(domz2(:,:,t+1))];
    spvmag(t+1) = norm(spv(:,t+1));
    
    %------------------ renormalization and LLE ----------------------
    if rem(t,normtime) == 0
        spvstar = spv(:,t+1) - dot(spv(:,t+1), fpv(:,t+1)) .* fpv(:,t+1);
        lam2inst = [lam2inst, ...
            (1/tnorm) * log(abs(norm(spvstar)))];
        spv(:,t+1) = spvstar./norm(spvstar);
        dpsi2(:,:,t+1) = vecToLat(spv(1:N,t+1),Nx,Ny);
        domz2(:,:,t+1) = vecToLat(spv(N+1:2*N,t+1),Nx,Ny);
    end
    
    dpsi2(:,:,t+1) = vecToLat(spv(1:N,t+1),Nx,Ny);
    domz2(:,:,t+1) = vecToLat(spv(N+1:2*N,t+1),Nx,Ny);

    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
    
end

lam2 = mean(lam2inst);

end