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

    % % ------------- predictor -------------
    dpsihguess = dpsihmat(:,:,k);
    domzhguess = domzhmat(:,:,k);
    for iters = 3
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

    %------ETD1 instead of predictor corrector
    % phi1L1 = (expm(p.L1*dt) - 1) ./ (p.L1+eps);
    % phi1L2 = (expm(p.L2*dt) - 1) ./ (p.L2+eps);
    % dpsihmat(:,:,k) = fft2(real(ifft2(dpsihmat(:,:,k)))) .* expm(p.L1*dt) ...
    %         + dt * phi1L1 .* N1hfunc(p,psihmat,dpsihmat(:,:,k),omzhmat,domzhmat(:,:,k));
    % domzhmat(:,:,k) = fft2(real(ifft2(domzhmat(:,:,k)))) .* expm(p.L2*dt) ...
    %         + dt * phi1L2 .* N2hfunc(p,psihmat,dpsihmat(:,:,k),omzhmat,domzhmat(:,:,k));

    %------symmetric usage is probably not necessary
    % dpsimat(:,:,k) = real(ifft2(dpsihmat(:,:,k),'symmetric'));
    % domzmat(:,:,k) = real(ifft2(domzhmat(:,:,k),'symmetric'));
    % dzetamat(:,:,k) = real(ifft2(dzetahmat(:,:,k),'symmetric'));
    dpsimat(:,:,k) = real(ifft2(dpsihmat(:,:,k)));
    domzmat(:,:,k) = real(ifft2(domzhmat(:,:,k)));
    dzetamat(:,:,k) = real(ifft2(dzetahmat(:,:,k)));
    
    dH(:,k) = [reshape(real(dpsimat(:,:,k)'),[],1);
                    reshape(real(domzmat(:,:,k)'),[],1)];
    % dH(:,k) = [latToVec(real(dpsimat(:,:,k)));
    %                 latToVec(real(domzmat(:,:,k)))];
    dHmag(k,1) = norm(dH(:,k));
    % dumat = real(ifft2(1i*Ky.*dzetahmat,'symmetric'));
    % dvmat = real(ifft2(-1i*Kx.*dzetahmat,'symmetric'));
end
end