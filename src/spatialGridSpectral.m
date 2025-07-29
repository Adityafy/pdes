function smesh = spatialGridSpectral(Lx,Nx)
% smesh = spatialGridSpectral(Lx,Nx)
% Generates the wavenumber grid for spectral and pseudospectral
% calculations. kx and ky are vectors, Kx and Ky makeup the wavenumber
% grid.
Ly = Lx;
Ny = Nx;
kx = (2*pi/Lx)*(-Nx/2:Nx/2-1)';
ky = (2*pi/Ly)*(-Ny/2:Ny/2-1)';
kx = fftshift(kx);
ky = fftshift(ky);
[Kx,Ky] = meshgrid(kx,ky);

smesh = struct('kx',kx,'ky',ky,'Kx',Kx,'Ky',Ky);
end