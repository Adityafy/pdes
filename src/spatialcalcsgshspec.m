function sp = spatialcalcsgshspec(Lx,cp)

Ly = Lx;
Nx = Lx;
Ny = Nx;
N = Nx*Ny;

% ---------  spatial grid in real domain -----------
x2 = (2*pi/16)*linspace(-Lx/2,Lx/2,Nx+1)';
x = x2(1:Nx);
y2 = (2*pi/16)*linspace(-Ly/2,Ly/2,Ny+1)';
y = y2(1:Ny);
[X,Y] = meshgrid(x,y);

% Linear derivative operator for Fourier domain
% for both the equations.

% wavenumber grid
kx = (2*pi/Nx)*(-Nx/2:Nx/2-1)';
ky = (2*pi/Ny)*(-Ny/2:Ny/2-1)';
% kx(Nx/2+1) = sqrt(eps);
% ky(Ny/2+1) = sqrt(eps);
kx = fftshift(kx);
ky = fftshift(ky);
[Kx,Ky] = meshgrid(kx,ky);

epsilon = cp.epsilon;
sigma = cp.sig;
c = sqrt(cp.csq);

L1 = epsilon - 1 - Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 ;
L2 = - sigma * (Kx.^2 + Ky.^2) - sigma * c^2;
Kdiff = Kx.^2 + Ky.^2;
firstrowKdiff = latToVec(Kdiff)';
KdiffMat = toeplitz(firstrowKdiff);
qx = latToVec(Kx);
qy = latToVec(Ky);

sp = struct('Lx',Lx,'Ly',Ly,'Nx',Nx,'Ny',Ny,'N',N,'x',x,'y',y,'X',X,'Y',Y, ...
    'kx',kx,'Kx',Kx,'Ky',Ky,'L1',L1,'L2',L2,'Kdiff',Kdiff,'qx',qx,'qy',qy);
end