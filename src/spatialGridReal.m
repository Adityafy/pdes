function rmesh = spatialGridReal(lam_0,Gamma)
% rmesh = spatialGridReal(lam_0,wlmult)
% lam_0 is the critical wavenumber
% wlmult is what I call the wavelength multiplier
% wlmult is something I came up with which is set-up in a way so that it
% eventually determines the Lx (length of the domain) to be multiple of
% wavelengths (wlmult multiples) and it determines delta x to be 2*pi/8 by
% setting Nx to be 8 times the wlmult.
% - Aditya

wlmult = Gamma/2;
Lx = lam_0 * wlmult;
Ly = Lx;
Nx = 8 * wlmult;
Ny = Nx;
N = Nx*Ny;

% ---------  spatial grid in real domain -----------
x2 = linspace(-Lx/2,Lx/2,Nx+1)';
x = x2(1:Nx);
y2 = linspace(-Ly/2,Ly/2,Ny+1)';
y = y2(1:Ny);
dx = x(2) - x(1);
dy = y(2) - y(1);

[X,Y] = meshgrid(x,y);

rmesh = struct('Lx',Lx,'Ly',Ly,'Nx',Nx,'Ny',Ny,'N',N,'x',x,'y',y, ...
                'dx',dx,'dy',dy,'X',X,'Y',Y);
end