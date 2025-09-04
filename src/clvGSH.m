function [wf,cle] = clvGSH(p,psi,omz,zeta,gsVec,cmat)
% [wf,cle] = clv(paramVec,dynamicsVec,v,cmat,couplingMat,jacobianFunc)

tu = size(cmat,3);
nv = p.ts.nv;
N = p.rmesh.N;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;
wi = zeros(size(gsVec)); % w (clv) initial
% wmag(:,1) = zeros(Nx,1); % w magnitude
cle(:,1) = zeros(nv,1);
nnorm = p.ts.nnorm;

fprintf("\nCalculating CLVs...\n");
for n = 1:tu-1
    for j = 1:nv
        for k = 1:j
            wi(:,j,n) = wi(:,j,n) + cmat(k,j,n).*gsVec(:,k,n);
        end
    end
    if rem(n,p.ts.nnorm) == 0
        for j = 1:nv
            wi(:,j,n) = wi(:,j,n)./norm(wi(:,j,n));
        end
    end
    if rem(n, tu/10) == 0
        fprintf("This is time step %g and ", n);
        toc
    end
end
figure;
wf = zeros(2*N, nv, tu); % w final
for n = 1:tu-1
    % zetamat = zetaGSHspectral(p,omzmat);
    % psi, omz,  pseudospectral, etd
    % [psimat,omzmat,zetamat] = ...
    %     advGSHstepPSETD(p,psimat,omzmat,zetamat,...
    %     etd.expL1dtmat,etd.expL2dtmat,@N1hat,@N2hat);

    %================== TS ==================
    for k = 1:nv
        dpsimat(:,:,k) = reshape(wi(1:N,k,n),Nx,Ny)';
        domzmat(:,:,k) = reshape(wi(N+1:2*N,k,n),Nx,Ny)';
        dzetamat(:,:,k) = zetaGSHspectral(p,domzmat(:,:,k));
    end
    
    % zeta(:,:,n) = zetaGSHspectral(p,omz(:,:,n));
    [dpsimat,domzmat,dzetamat,wfn,wfmagn] = ...
        advGSHTSstepPSETD(p,psi(:,:,n),omz(:,:,n),zeta(:,:,n),...
        dpsimat,domzmat,dzetamat,...
        p.etd.expL1dtmat,p.etd.expL2dtmat,@R1ts,@R2ts);

    % --------- renormalization and CLE ----------
    % if rem(n,nnorm) == 0
    %     idx = round(n/nnorm);
        for k = 1:nv
            wfmagn(k,1) = norm(wfn(:,k));
            cle(k,n) = (1/p.ts.tN) * log(abs(wfmagn(k,1)));
            wfn(:,k) = wfn(:,k)./norm(wfn(:,k));
        end
    % end
    wf(:,:,n) = wfn;
    % wfmagn(k,n) = norm(wf(:,k,n));
    imagesc(reshape(wf(1:p.rmesh.N,1,n),p.rmesh.Nx,p.rmesh.Ny)'); 
    colorbar; colormap jet; axis square; clim([-0.05 0.05]);
    drawnow;
end
toc
end
