function [psimat,omzmat,zetamat,umat,vmat] = ...
                    advGSHstepPSETD(p,psimat,omzmat,zetamat,...
                                    expL1hmat,expL2hmat,N1hfunc,N2hfunc)
    dt = p.sim.dt;
    Kx = p.smesh.Kx;
    Ky = p.smesh.Ky;

    psihmat = fft2(psimat);
    omzhmat = fft2(omzmat);
    zetahmat = fft2(zetamat);
    % ------------- predictor -------------
    psihguess = psihmat;
    omzhguess = omzhmat;
    for iters = 2
        psihpred = fft2(real(ifft2(psihmat))) .* expL1hmat ...
                    + 0.5 * dt * N1hfunc(p,psihguess,zetahmat) ...
                    + 0.5 * dt * expL1hmat .* N1hfunc(p,psihmat,zetahmat);
        omzhpred = fft2(real(ifft2(omzhmat))).*expL2hmat ...
                    + 0.5 * dt * N2hfunc(p,psihguess) ...
                    + 0.5 * dt * expL2hmat.*N2hfunc(p,psihmat);
        psihguess = psihpred;
        omzhguess = omzhpred;
    end
    % -------------corrector -------------
    psihmat = fft2(real(ifft2(psihmat))) .* expL1hmat ...
                    + 0.5 * dt * N1hfunc(p,psihpred,zetahmat) ...
                    + 0.5 * dt * expL1hmat .* N1hfunc(p,psihmat,zetahmat);

    omzhmat = fft2(real(ifft2(omzhmat))) .* expL2hmat ...
                    + 0.5 * dt * N2hfunc(p,psihpred) ...
                    + 0.5 * dt * expL2hmat.*N2hfunc(p,psihmat);

    psimat = real(ifft2(psihmat,'symmetric'));
    omzmat = real(ifft2(omzhmat,'symmetric'));
    zetamat = real(ifft2(zetahmat,'symmetric'));
    
    umat = real(ifft2(1i*Ky.*zetahmat,'symmetric'));
    vmat = real(ifft2(-1i*Kx.*zetahmat,'symmetric'));
end