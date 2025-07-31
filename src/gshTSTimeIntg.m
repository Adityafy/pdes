function [psi,omz,zeta,u,v,pertvecs,dHmag,laminst,lamgs] = gshTSTimeIntg(p,dynICAddress)

% dt = p.ts.dt;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;
N = p.rmesh.N;
nv = p.ts.nv;
tN = p.ts.tN;
nmax = p.sim.nmax;
interv = p.sim.interv;
nnorm = p.ts.nnorm;
if p.sim.runtype == 0
        [matdivpsi, matdivomz] = impMatGSH(p,p.sim.dt);
end
[psimat,psi,omzmat,omz,zetamat,zeta] = dynICafterTr(dynICAddress);
[dpsimat, dpsihmat, domzmat, domzhmat, ...
                    dzetamat, dH, laminst, dHmag] = tsics(p);

fprintf('Initial ||dH(:,1)|| = %.4e\n', norm(dH(:,1)));
% figure;
% imagesc(dpsimat(:,:,1)); title('Initial \delta\psi_n^{(1)}'); 
% colormap jet; colorbar; axis square; drawnow;
% 
% figure;
% [psimat,psi,omzmat,omz,zetamat,zeta,u,v] = ICfunc(p,p.sim.dynICtype);
for n = 1:nmax
    %================== saving dynamics at intervals ==================
    if p.sim.runtype == 0
        [psimat,omzmat,zetamat,umat,vmat] = advGSHstepFDSI(psimat, ...
            omzmat, zetamat, matdivpsi, matdivomz, p);
        [psimat, omzmat, dH1, dHmag, ~] = advGSHdH1stepFDSI( ...
            psimat, omzmat, zetamat, dpsimat(:,:,1), domzmat(:,:,1), dzetamat(:,:,1), ...
            matdivpsi, matdivomz, n, dHmag, 0, 0, p);
        %------------------ renormalization and LLE ----------------------
        if rem(n,nnorm) == 0
            lam1inst = [lam1inst, (1/tN) * log(abs(norm(dH1)))];
            dH1 = dH1./norm(dH1);
            dpsi1mat = reshape(dH1(1:N),Nx,Ny)';
            domz1mat = reshape(dH1(N+1:2*N),Nx,Ny)';
        end
    
        dpsi1mat = reshape(dH1(1:N),Nx,Ny)';
        domz1mat = reshape(dH1(N+1:2*N),Nx,Ny)';

    elseif p.sim.runtype == 1
        % zeta,  pseudospectral
        zetamat = zetaGSHspectral(p,omzmat);
        % psi, omz,  pseudospectral, etd
        [psimat,omzmat,zetamat,umat,vmat] = ...
            advGSHstepPSETD(p,psimat,omzmat,zetamat,...
            p.etd.expL1dtmat,p.etd.expL2dtmat,@N1hat,@N2hat);

        %================== TS ==================
        [dpsimat,domzmat,dzetamat,dH,dHmagn] = ...
            advGSHTSstepPSETD(p,psimat,omzmat,zetamat,...
            dpsimat,domzmat,dzetamat,...
            p.etd.expL1dtmat,p.etd.expL2dtmat,@R1ts,@R2ts);

        % if rem(n/nnorm) == 0
        %     [dpsimat, domzmat, dpsihmat, domzhmat, laminst] = ...
        %         gsStep(p, n, dpsimat, domzmat, dpsihmat, domzhmat, dH, dHmagn);
        % end
        % --------- renormalization and LLE ----------
        if rem(n,nnorm) == 0
            %----just using one vector
            % dH(:,1) = dH(:,1) / norm(dH(:,1));
            % laminst(1,round(n/nnorm)) = (1/tN) * log(norm(dH(:,1)));
            
            %----using multiple vectors
            [Q,R]= qr(dH,'econ','vector');
            dH = Q(:,1:nv);
            dH1(:,round(n/nnorm)) = dH(:,1);
            for k = 1:nv
                laminst(k,round(n/nnorm)) = (1/tN)*log(abs(R(k,k)));
            end

            for k = 1:nv
                dpsimat(:,:,k) = reshape(dH(1:N,k),Nx,Ny)';
                domzmat(:,:,k) = reshape(dH(N+1:2*N,k),Nx,Ny)';
                % dpsimat(:,:,k) = vecToLat(dH(1:N,k),Nx,Ny);
                % domzmat(:,:,k) = vecToLat(dH(N+1:2*N,k),Nx,Ny);
                dpsihmat(:,:,k) = fft2(dpsimat(:,:,k));
                domzhmat(:,:,k) = fft2(domzmat(:,:,k));
                dHmagn(k,1) = norm(dH(:,k));
            end
        end

    end
    %================== saving dynamics at intervals ==================
    dHmag(:,n) = dHmagn;
    % laminst(:,n) = laminstn;
    if rem(n,nmax/interv) == 0
        psi(:,:,round(n*interv/nmax)) = psimat;
        omz(:,:,round(n*interv/nmax)) = omzmat;
        zeta(:,:,round(n*interv/nmax)) = zetamat;
        u(:,:,round(n*interv/nmax)) = umat;
        v(:,:,round(n*interv/nmax)) = vmat;
        pertvecs(:,:,round(n*interv/nmax)) = dH;
        % dpsi1(:,:,round(n*interv/nmax)) = dpsi1mat;
        % domz1(:,:,round(n*interv/nmax)) = domz1mat;
    end
    if rem(n,round(nmax/10)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax);
        toc;
    end
    if isnan(abs(sum(sum(psimat))))
        error('blow up occured in psi..');
    end
    if abs(sum(sum(dpsimat(:,:,1)))) < 10^-5
        fprintf('time step : %g, ', n);
        error('Perturbation vector elements went to zero');
    end
    
    if p.sim.makeLiveFig == 1
        subplot(1,2,1);
        imagesc(psimat); colorbar; colormap jet; clim([-1 1]); axis square;
        % imagesc([psimat psimat]); colorbar; colormap jet; clim([-0.8 0.8]);
        subplot(1,2,2);
        imagesc(dpsimat(:,:,1));
        % imagesc(dpsi1mat);
        colorbar; axis square;  colormap jet; clim([-0.1 0.1]);
        % clim([-1 1]);
        % % subplot(1,3,3); imagesc(zetatrmat); colorbar; axis square;
        drawnow;
    end
end

for k = 1:nv
    lamgs(k,1) = (1/(p.sim.tu/p.ts.tN))*sum(laminst(k,:));
end

end