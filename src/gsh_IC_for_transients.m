function [psimat,psi,omzmat,omz,zetamat,zeta,u,v] = gsh_IC_for_transients(p,ictype)
% [psimat,psi,omzmat,omz,zetamat,zeta,u,v] = gsh_IC_for_transients(p,ictype)
% INITIAL CONDITIONS
% make sure zeta (zetatrmat) initialization is not too high (max < 0.1) or
% the code breaks.
% ictype = 0 for periodic, 1 for random

seed = p.sim.seed;
X = p.rmesh.X;
Y = p.rmesh.Y;
Kx = p.smesh.Kx;
Ky = p.smesh.Ky;

if ictype == 0
    % -----periodic ICs-----
    %%% note: the angular frequency should be an integral multiple of 2*pi/L
    psimat = 0.01*(cos(6*pi*X/Lx)+sin(6*pi*Y/Ly));
    zetamat = 0.01*(cos(2*pi*X/Lx)+sin(2*pi*Y/Ly));
elseif ictype == 1
    % -----random ICs-----
    rng(seed,"twister");
    psimat = -0.1 + (0.1+0.1) * rand(length(X), length(Y));
    rng(seed + 1,"twister");
    zetamat = -0.1 + (0.1+0.1)*rand(length(X), length(Y));
    rng(seed,"twister");
end


%%%======= Initializing variables for time marching ========
psi(:,:,1) = psimat; % 3D mat that's going to be appended and saved
% psitrvec = reshape(psitr',[],1); % changing to vector
% psitrvec = fftshift(psitrvec);
psihmat = fft2(psimat); % assigning psihat for fft of psi

zeta = zetamat;
umat = zeros(length(X), length(Y),1);
vmat = zeros(length(X), length(Y),1);
u = zeros(length(X), length(Y),1);
v = zeros(length(X), length(Y),1);
% zetatrvec = reshape(zetatrmat',[],1);
% zetatrvec = fftshift(zetatrvec);
zetahmat = fft2(zeta);

% omztr = -fdlaplacian(zetatr, Nx,Ny,dx); % if fd laplacian is to be used
omzmat = (Kx.^2+Ky.^2).*(zetamat);
omz = omzmat;
omzhmat = fft2(omzmat);
end