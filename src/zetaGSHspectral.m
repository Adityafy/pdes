function zetamat = zetaGSHspectral(p,omzmat)
% zetahmat = zetaPS(omzmat,Kdiff)
% Calculates the Poisson equation in GSH using spectral method.
Kdiff = p.Kdiff;
Kdiff(1,1) = 1;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;
zetamat2h = fft2(omzmat)./(Kdiff);
zetamat2h(1,1) = 0+1i*0;
% zetamat2h(1,Ny) = 0+1i*0;
% zetamat2h(Nx,1) = 0+1i*0;
% zetamat2h(Nx,Ny) = 0+1i*0;
zetahmat = fft2(real(ifft2(zetamat2h)));
zetamat = real(ifft2(zetahmat));
% imagesc(zetamat); colorbar; axis square;
end