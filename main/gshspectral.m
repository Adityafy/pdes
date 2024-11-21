%% Script for GSH with Spectral
% Using ETD1 explicit (predictor corrector)

clear all;
close all;

addpath('../src/');
mfa = '../../pdesDataDump/media/'; % media folder address
dfa = '../../pdesDataDump/data/'; % saved data folder address

%% Parameter struct

% key parameters
eps = 0.7;
sig = 1;
csq = 0.2;
c = sqrt(csq);
gm = 50;
lam_0 = 2*pi; % critical roll wavelength
dt_tr = 0.2;
dt_fpv = 0.001;
gamma = 16;
trtu = 250; % transient time units
totu = 0; % total time units
seed = 4;
tN = 2; % renormalization time units
spatial_res = 16;
dx = lam_0/spatial_res; % node spacing
dy = dx;
Nx = gamma*4;
Ny = Nx;
N = Nx*Ny;

totmax = totu/dt_fpv;
nmax = trtu/dt_tr;

% spatial grid in real domain
x = spatial_res*2*pi*(1:Nx)'/Nx;
y = spatial_res*2*pi*(1:Ny)'/Ny;

[X,Y] = meshgrid(x,y);

%% Linear derivative operator for Fourier domain
% for both the equations.

% wavenumber grid
kx = (1/spatial_res)*(-Nx/2:Nx/2-1)';
ky = (1/spatial_res)*(-Ny/2:Ny/2-1)';
kx(Nx/2+1) = 10^(-5);
ky(Ny/2+1) = 10^(-5);
kx = fftshift(kx);
ky = fftshift(ky);
[Kx,Ky] = meshgrid(kx,ky);

L1 = eps - 1 - Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 ;
L2 = - sig * (Kx.^2 + Ky.^2) - sig * c^2;

%% Intervals
TrIntervals = floor(trtu);
ToIntervals = floor(totu);
if totu < 1
    ToIntervals = 10;
end

%% INITIAL CONDITIONS

% -----periodic ICs-----
% u = cos(x/spatial_res).*(1+sin(x/spatial_res)); 

% -----random ICs-----
rng(seed,"twister");
psi(:,:,1) = -0.1 + (0.1+0.1) * rand(length(X), length(Y));
psimat = -0.1 + (0.1+0.1) * rand(length(X), length(Y));
% changing to vector
psivec = latToVec(psi);
% assigning psihat for fft of psi
psivec = fftshift(psivec);
psih = fft(psivec); % psihat

%-------------- random IC for zeta and omega -----------------
zetamat = zeros(length(X), length(Y),1); % scalar field
umat = zeros(length(X), length(Y),1);
vmat = zeros(length(X), length(Y),1);
u = zeros(length(X), length(Y),1);
v = zeros(length(X), length(Y),1);
zetavec = latToVec(zetamat);
zetavec = fftshift(zetavec);
zetah = fft(zetavec);

rng(seed + 1,"twister");
omz = -1 + (0.1+0.1)*rand(length(X), length(Y));
omzmat = omz;
omzvec = latToVec(omz);
omzvec = fftshift(omzvec);
omzh = fft(omzvec);


%% ETD explicit
% v = fft(u(:,1));
% tic
% psi = etd1exp_sh2D(psi,v,dt_tr,nmax,L,ToIntervals,Nx,Ny);
qx = latToVec(Kx);
qy = latToVec(Ky);
expL1hvec = exp(latToVec(L1*dt_tr));
expL2hvec = exp(latToVec(L2*dt_tr));

tic;
for n = 1:nmax
    %------------------ zeta, iterative ------------------------
    [zetamat,difference] = iterZetaOmz(omzmat,zetamat,X(1,2)-X(1,1),Nx,Ny);
    zetavec = latToVec(zetamat);
    zetah = fft(zetavec);
    % zetamath = fft2(zetamat);
    % zetah = latToVec(zetamath);

    % zetah = fft(real( ifft(omzh./(qx.^2+qy.^2)) ));
    % zetamath = fft2(omzmat)./(Kx.^2+Ky.^2);
    % zetamat = ifft2(zetamath);
    % zetah = fft(real(latToVec(zetamat)));

    % -------------predictor-------------
    % psihguess = fft(latToVec(0.1*rand(length(X),length(Y)))); 
                      % (initial guess for predictor step
                        % This guess is going to updated to psipred so
                        % initializing that already)
    psihguess = psih;
    % omzhguess = fft(latToVec(rand(length(X),length(Y))));
    omzhguess = omzh;
    for iters = 2
        psihpred = fft(real(ifft(psih))) .* expL1hvec ...
                    + 0.5 * dt_tr * N1hat(psihguess,zetah,qx,qy) ...
                    + 0.5 * dt_tr * expL1hvec .* N1hat(psih,zetah,qx,qy);
        omzhpred = fft(real(ifft(omzh))).*expL2hvec ...
                    + 0.5 * dt_tr * N2hat(gm,psihguess,qx,qy) ...
                    + 0.5 * dt_tr * expL2hvec.* N2hat(gm,psih,qx,qy);
        psihguess = psihpred;
    end
    % -------------corrector, psi -------------
    psih = fft(real(ifft(psih))) .* expL1hvec ...
                    + 0.5 * dt_tr * N1hat(psihpred,zetah,qx,qy) ...
                    + 0.5 * dt_tr * expL1hvec .* N1hat(psih,zetah,qx,qy);

    % -------------corrector, omz -------------
    omzh = fft(real(ifft(omzh))) .* expL2hvec ...
                    + 0.5 * dt_tr * N2hat(gm,psihpred,qx,qy) ...
                    + 0.5 * dt_tr * expL2hvec.* N2hat(gm,psih,qx,qy);
    

    psimat = vecToLat(real(ifft(psih)),Nx,Ny);
    omzmat = vecToLat(real(ifft(omzh)),Nx,Ny);
    zetamat = vecToLat(real(ifft(zetah)),Nx,Ny);
    
    umat = vecToLat(real(ifft(1i*qy.*zetah)),Nx,Ny);
    vmat = vecToLat(real(ifft(-1i*qx.*zetah)),Nx,Ny);

    % psi(:,:,n+1) = psimat;
    % omz(:,:,n+1) = omzmat;
    %------------- saving dynamics at intervals ---------------
    if rem(n,nmax/TrIntervals) == 0
        psi(:,:,round(n*TrIntervals/nmax)) = psimat;
        omz(:,:,round(n*TrIntervals/nmax)) = omzmat;
        zeta(:,:,round(n*TrIntervals/nmax)) = zetamat;
        u(:,:,round(n*TrIntervals/nmax)) = umat;
        v(:,:,round(n*TrIntervals/nmax)) = vmat;
    end
    if rem(n,round(nmax/10)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax);
        toc;
    end
    % figure;
    % plot(difference,'-o');
    % 
    % figure;
    % quiver(X,Y,vmat,umat,3);
    % 
    % figure;
    % contourf(omzmat,'LevelStep',0.05,'EdgeColor','none'); colorbar;
    % close all;

end

%%
% figure; contourf(psi(:,:,4),'LevelStep',0.01,'EdgeColor','none');
% colormap jet; colorbar;
% set(gca,'YDir','normal'); axis square;

%%
figure; hold on;
attime = 1;
contourf(psi(:,:,attime),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
[X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,vmat(:,:,attime),umat(:,:,attime),3,'black');
hold off;

%%
figure;
subplot(1,2,1);
hold on;
attime = n;
contourf(psi(:,:,end),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
% [X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,vmat(:,:,end),umat(:,:,end),3,'black');
hold off;
subplot(1,2,2);
contourf(omz(:,:,end),'LevelStep',0.01,'EdgeColor','none');
colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;


%% functions

function N1h = N1hat(psih,zetah,qx,qy)
% function for the nonlinear part of equation 1
N1h = -fft(real((ifft(psih)).^3)) ...
        - fft(real( ifft(1i*qy.*zetah) .* ifft(1i*qx.*psih) )) ...
            + fft(real( ifft(1i*qx.*zetah) .* ifft(1i*qy.*psih) ));
end

function N2h = N2hat(gm,psih,qx,qy)
% function for the nonlinear part of equation 2

% derivatives in the fourier domain
psi_y = ifft(1i*qy.*psih);
psi_xxx = ifft(-1i*qx.^3.*psih);
psi_xyy = ifft(-1i*qx.*qy.^2.*psih);
psi_x = ifft(1i*qx.*psih);
psi_yxx = ifft(-1i*qy.*qx.^2.*psih);
psi_yyy = ifft(-1i*qy.^3.*psih);

% fft of the nonlinear part
N2h = fft(real(-gm * psi_y.*(psi_xxx+psi_xyy))) ...
     + fft(real(gm * psi_x.*(psi_yxx+psi_yyy)));
end

function [zeta,difference] = iterZetaOmz(omzmat, zetamat, dx, Nx, Ny)
% zeta = iterativeZetaOmz(omzmat, zetamat, p)
% An iterative poisson solver for
% Omega_z = - Laplacian(Zeta)
% with periodic boundary conditions.
% mind the negative sign if used for other kinds of Poisson solvers
% Number of iterations is fixed for now (= 800)
    % delta = p.dx;

    I = 1:Nx;
    J = 1:Ny;
    Ip1 = circshift(I,-1);
    Im1 = circshift(I,1);
    Jp1 = circshift(J,-1);
    Jm1 = circshift(J,1);

    zeta = zetamat;
    iterations = 400;
    difference = zeros(1,iterations);

    for k = 1:iterations
        for i = 1:length(I)
            for j = 1:length(J)
                zeta(i,j) = ((dx^2)/4) * omzmat(I(i),J(j)) ...
                    + (1/4) * ( zetamat(Im1(i),J(j)) + zetamat(Ip1(i),J(j)) ...
                    + zetamat(I(i),Jm1(j)) + zetamat(I(i),Jp1(j)) );
            end
        end

        lapl = zeros(size(zeta));
        for i = 1:length(I)
            for j = 1:length(J)
                lapl(i,j) = (1/(dx^2)) * ...
                    ( zeta(Im1(i),J(j)) + zeta(Ip1(i),J(j)) ...
                    - 4 * zeta(I(i),J(j)) + ...
                    + zeta(I(i),Jm1(j)) + zeta(I(i),Jp1(j)) );
            end
        end
        b = lapl - (-omzmat);
        difference(k) = max(max(abs(b)));
        zetamat = zeta;
    end
end
