function rhshat = R1ts(p,psihmat,dpsihmat,omzhmat,domzhmat)

Kx = p.smesh.Kx;
Ky = p.smesh.Ky;
Kdiff = p.Kdiff;
frgamma = p.sim.frgamma;
frgamma = 1;

omzmat = real(ifft2(omzhmat));
domzmat = real(ifft2(domzhmat));

zetamat = zetaGSHspectral(p,omzmat);
zetahmat = fft2(zetamat);

dzetamat = zetaGSHspectral(p,domzmat);
dzetahmat = fft2(dzetamat);

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