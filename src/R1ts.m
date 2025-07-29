function rhshat = R1ts(p,psihmat,dpsihmat,omzhmat,domzhmat)

Kx = p.smesh.Kx;
Ky = p.smesh.Ky;
Kdiff = p.Kdiff;
frgamma = p.sim.frgamma;
frgamma =1;

% psihmat = fft2(psimat);
% dpsihmat = fft2(dpsimat);
% omzhmat = fft2(omzmat);
% domzhmat = fft2(domzmat);
% zetahmat = omzhmat./Kdiff;
% zetahmat(1,1) = 0;
% dzetahmat = domzhmat./Kdiff;
% dzetahmat(1,1) = 0;

zetahmat = zetaGSHspectral(p,omzhmat);
dzetahmat = zetaGSHspectral(p,domzhmat);

term1 = -3*(ifft2(psihmat).^2) .* ifft2(dpsihmat);
unfiltered2 = - ifft2(1i*Ky.*zetahmat) .* ifft2(1i*Kx.*dpsihmat)  ...
    + ifft2(1i*Kx.*zetahmat) .* ifft2(1i*Ky.*dpsihmat) ...
    - ifft2(1i*Ky.*dzetahmat) .* ifft2(1i*Kx.*psihmat) ...
    + ifft2(1i*Kx.*dzetahmat) .* ifft2(1i*Ky.*psihmat) ;
filtered2 = imgaussfilt(real(unfiltered2),frgamma, ...
                        'FilterDomain','spatial','Padding','circular');
rhshat = fft2(real(term1) + filtered2);
% rhshat = fft2(real(term1) + unfiltered2);
end