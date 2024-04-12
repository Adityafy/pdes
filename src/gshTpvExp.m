function [dpsi3, domz3, dzeta3, tpv, tpvmag, lam3inst, lam3, res3] = ...
gshTpvExp(p,psi,omz,zeta,fpv,spv,evofunc1,evofunc2,fpvicfunc)
% [dpsi3, domz3, dzeta3, tpv, tpvmag, lam3inst, lam3, res3] = ...
%    gshTpvExp(p,psi,omz,zeta,fpv,spv,evofunc1,evofunc2,fpvicfunc)
% Evolution of the third perturbation vector (tpv) that uses an explicit
% method. 
Nx = p.Nx;
Ny = p.Ny;
dt = p.dt;
tmax = p.tmax;
N = Nx*Ny;

rng(p.seed+2,"twister");

[dpsi3,domz3,tpv] = fpvicfunc(p);

tpvmag = norm(tpv(:,1));

dzeta3 = zeros(Nx,Ny,tmax);

normtime = p.normtime;
tnorm = normtime*dt;
lam3inst = [];
res3 = zeros(1, tmax);

fprintf('\nEvolving the third perturbation vector...\n');

for t = 1:tmax
    
    %------------------ explicit TS linear ----------------------
    dzeta3(:,:,t) = iterativeZetaOmz(domz3(:,:,t),dzeta3(:,:,t),p);
    [dpsi3(:,:,t+1), res3(t)] = evofunc1(psi(:,:,t),zeta(:,:,t),dpsi3(:,:,t), ...
                        dzeta3(:,:,t),p,@ts1exp,@sfdd);
    domz3(:,:,t+1) = evofunc2(omz(:,:,t),psi(:,:,t),domz3(:,:,t), ...
                        dpsi3(:,:,t),p,@ts2exp,@sfdd);
    tpv(:,t+1) = [latToVec(dpsi3(:,:,t+1)); latToVec(domz3(:,:,t+1))];
    tpvmag(t+1) = norm(tpv(:,t+1));
    
    %------------------ renormalization and LLE ----------------------
    if rem(t,normtime) == 0
        tpvstar = tpv(:,t+1) - dot(tpv(:,t+1), spv(:,t+1)) .* spv(:,t+1) ...
                             - dot(tpv(:,t+1), fpv(:,t+1)) .* fpv(:,t+1);
        lam3inst = [lam3inst, ...
            (1/tnorm) * log(abs(norm(tpvstar)))];
        tpv(:,t+1) = tpvstar./norm(tpvstar);
        dpsi3(:,:,t+1) = vecToLat(tpv(1:N,t+1),Nx,Ny);
        domz3(:,:,t+1) = vecToLat(tpv(N+1:2*N,t+1),Nx,Ny);
    end
    
    dpsi3(:,:,t+1) = vecToLat(tpv(1:N,t+1),Nx,Ny);
    domz3(:,:,t+1) = vecToLat(tpv(N+1:2*N,t+1),Nx,Ny);

    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
    
end

lam3 = mean(lam3inst);

end