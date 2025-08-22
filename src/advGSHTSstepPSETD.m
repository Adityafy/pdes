function [dpsimat,domzmat,dzetamat,dH,dHmag] = ...
    advGSHTSstepPSETD(p,psimat,omzmat,zetamat,...
    dpsimat,domzmat,dzetamat,...
    expL1dtmat,expL2dtmat,N1hfunc,N2hfunc)

dt = p.sim.dt;
nv = p.ts.nv;
N = p.rmesh.N;
Kx = p.smesh.Kx;
Ky = p.smesh.Ky;

psihmat = fft2(psimat);
omzhmat = fft2(omzmat);

dpsihmat = zeros(size(dpsimat));
domzhmat = zeros(size(domzmat));
dzetahmat = zeros(size(dzetamat));
dH = zeros(2*N,nv);

for k = 1:nv
    dpsihmat(:,:,k) = fft2(dpsimat(:,:,k));
    domzhmat(:,:,k) = fft2(domzmat(:,:,k));
    dzetahmat(:,:,k) = fft2(dzetamat(:,:,k));

    % % ------------- predictor -------------
    dpsihguess = dpsihmat(:,:,k);
    domzhguess = domzhmat(:,:,k);
    for iters = 1
        dpsihpred = dpsihmat(:,:,k) .* expL1dtmat ...
            + 0.5 * dt * N1hfunc(p,psihmat,dpsihguess,omzhmat,domzhmat(:,:,k)) ...
            + 0.5 * dt * expL1dtmat .* N1hfunc(p,psihmat,dpsihguess,omzhmat,domzhmat(:,:,k));
        domzhpred = domzhmat(:,:,k) .* expL2dtmat ...
            + 0.5 * dt * N2hfunc(p,psihmat,dpsihguess,omzhmat,domzhmat(:,:,k)) ...
            + 0.5 * dt * expL2dtmat.*N2hfunc(p,psihmat,dpsihguess,omzhmat,domzhmat(:,:,k));
        dpsihguess = dpsihpred;
        domzhguess = domzhpred;
    end
    % -------------corrector -------------
    dpsihmat(:,:,k) = dpsihmat(:,:,k) .* expL1dtmat ...
        + 0.5 * dt * N1hfunc(p,psihmat,dpsihpred,omzhmat,domzhmat(:,:,k)) ...
        + 0.5 * dt * expL1dtmat .* N1hfunc(p,psihmat,dpsihpred,omzhmat,domzhmat(:,:,k));

    domzhmat(:,:,k) = domzhmat(:,:,k) .* expL2dtmat ...
        + 0.5 * dt * N2hfunc(p,psihmat,dpsihpred,omzhmat,domzhmat(:,:,k)) ...
        + 0.5 * dt * expL2dtmat.*N2hfunc(p,psihmat,dpsihpred,omzhmat,domzhmat(:,:,k));

    % going back to real space
    dpsimat(:,:,k) = real(ifft2(dpsihmat(:,:,k)));
    domzmat(:,:,k) = real(ifft2(domzhmat(:,:,k)));
    dzetamat(:,:,k) = real(ifft2(dzetahmat(:,:,k)));
    
    dH(:,k) = [reshape(real(dpsimat(:,:,k)'),[],1);
                    reshape(real(domzmat(:,:,k)'),[],1)];

    dHmag(k,1) = norm(dH(:,k));

end
end