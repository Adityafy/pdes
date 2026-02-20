function [umat, vmat] = uvZetaSpectral(p,zetamat)
Kx = p.smesh.Kx;
Ky = p.smesh.Ky;
zetahmat = fft2(zetamat);
umat = real(ifft2(1i*Ky.*zetahmat,'symmetric'));
vmat = real(ifft2(-1i*Kx.*zetahmat,'symmetric'));
end