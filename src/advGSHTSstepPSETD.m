function [dpsimat,domzmat,dzetamat,dH,dHmag] = ...
    advGSHTSstepPSETD(p,psimat,omzmat,zetamat,...
    dpsimat,domzmat,dzetamat,...
    expL1dtmat,expL2dtmat,N1hfunc,N2hfunc)

dt = p.sim.dt;
nv = p.ts.nv;
Kx = p.smesh.Kx;
Ky = p.smesh.Ky;

psihmat = fft2(psimat);
omzhmat = fft2(omzmat);
% zetahmat = fft2(zetamat);

for k = 1:nv
    dpsihmat(:,:,k) = fft2(dpsimat(:,:,k));
    domzhmat(:,:,k) = fft2(domzmat(:,:,k));
    dzetahmat(:,:,k) = fft2(dzetamat(:,:,k));
    % ------------- predictor -------------
    dpsihguess = dpsihmat(:,:,k);
    domzhguess = domzhmat(:,:,k);
    for iters = 2
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

    dpsimat(:,:,k) = real(ifft2(dpsihmat(:,:,k),'symmetric'));
    domzmat(:,:,k) = real(ifft2(domzhmat(:,:,k),'symmetric'));
    dzetamat(:,:,k) = real(ifft2(dzetahmat(:,:,k),'symmetric'));
    
    dH(:,k) = [reshape(real(dpsimat(:,:,k)'),[],1);
                    reshape(real(domzmat(:,:,k)'),[],1)];
    dHmag(k,1) = norm(dH(:,k));
    % dumat = real(ifft2(1i*Ky.*dzetahmat,'symmetric'));
    % dvmat = real(ifft2(-1i*Kx.*dzetahmat,'symmetric'));
end
end