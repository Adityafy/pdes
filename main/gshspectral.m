%% Script for GSH with Spectral
% Using ETD1 explicit (predictor corrector)

clear all;
% close all;

addpath('../src/');
mfa = '../../pdesDataDump/media/'; % media folder address
dfa = '../../pdesDataDump/data/'; % saved data folder address

%% Parameter struct

% key parameters
epsilon = 0.7;
sig = 1;
csq = 0.2;
c = sqrt(csq);
gm = 50;
lam_0 = 2*pi; % critical roll wavelength
dt_tr = 0.1;
dt_fpv = 0.001;
% gamma = 15;
trtu = 100; % transient time units
totu = 0; % total time units
seed = 4;
tN = 2; % renormalization time units
% spatial_res = 8;
% dx = lam_0/spatial_res; % node spacing
% dy = dx;
% Nx = gamma*(spatial_res/4);
Lx = 64;
Ly = Lx;
Nx = Lx;
Ny = Nx;
N = Nx*Ny;

totmax = totu/dt_fpv;
nmax = trtu/dt_tr;

% ---------  spatial grid in real domain -----------
x2 = (2*pi/16)*linspace(-Lx/2,Lx/2,Nx+1)';
x = x2(1:Nx);
y2 = (2*pi/16)*linspace(-Ly/2,Ly/2,Ny+1)';
y = y2(1:Ny);
dx = x(2) - x(1);
dy = y(2) - y(1);
[X,Y] = meshgrid(x,y);

%% Linear derivative operator for Fourier domain
% for both the equations.

% wavenumber grid
kx = (2*pi/Nx)*(-Nx/2:Nx/2-1)';
ky = (2*pi/Ny)*(-Ny/2:Ny/2-1)';
% kx(Nx/2+1) = 1e-5;
% ky(Ny/2+1) = 1e-5;
kx = fftshift(kx);
ky = fftshift(ky);
[Kx,Ky] = meshgrid(kx,ky);

L1 = epsilon - 1 - Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 ;
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
% zetamat = cos(X/16).*(1+sin(Y/16));
zeta = zeros(size(zetamat));
umat = zeros(length(X), length(Y),1);
vmat = zeros(length(X), length(Y),1);
u = zeros(length(X), length(Y),1);
v = zeros(length(X), length(Y),1);
zetavec = latToVec(zetamat);
zetavec = fftshift(zetavec);
zetah = fft(zetavec);

rng(seed + 1,"twister");
omz = -0.1 + (0.1+0.1)*rand(length(X), length(Y));
omzmat = omz;
omzvec = latToVec(omz);
omzvec = fftshift(omzvec);
omzh = fft(omzvec);

%% coefficient matrix for poisson
firstrow = zeros(1,Nx*Ny);
firstrow(1) = 4/dx^2;
firstrow(2) = -1/dx^2;
firstrow(Nx) = -1/dx^2;
firstrow(Ny+1) = -1/dx^2;
firstrow(Nx*Ny-Nx+1) = -1/dx^2;
DiffMatZetaOmz = sparse(toeplitz(firstrow));

%% ETD explicit
Kdiff = Kx.^2 + Ky.^2;
firstrowKdiff = latToVec(Kdiff)';
KdiffMat = toeplitz(firstrowKdiff);
qx = latToVec(Kx);
qy = latToVec(Ky);
expL1hvec = exp(latToVec(L1*dt_tr));
expL2hvec = exp(latToVec(L2*dt_tr));


tic;

for n = 1:nmax
    
    % %------------------ zeta, pseudospectral ------------------------
    zetamat2h = fft2(omzmat)./(Kdiff+eps);
    zetah = fft(real(latToVec(ifft2(zetamat2h))));
    zetamat = ifft2(zetamat2h);
    poisson_test = isequal(omzmat,ifft2(Kdiff*zetamat2h));
    % zetah = fft(real(latToVec(ifft2(zeta2h))));

    % -------------predictor-------------
    % psihguess = fft(latToVec(0.1*rand(length(X),length(Y)))); 
    psihguess = psih;   % (initial guess for predictor step
                        % This guess is going to updated to psipred so
                        % initializing that already)
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
    if isnan(abs(sum(sum(psimat))))
        error('blow up occured in psi..');
    end

end

%%
% figure; contourf(psi(:,:,4),'LevelStep',0.01,'EdgeColor','none');
% colormap jet; colorbar;
% set(gca,'YDir','normal'); axis square;

%%
figure; hold on;
attime = trtu;
contourf(X,Y,psi(:,:,attime),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
% [X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,v(:,:,attime),u(:,:,attime),3,'black');
box on;
xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
hold off;

%%
figure;
subplot(1,2,1);
hold on;
attime = trtu;
contourf(X,Y,psi(:,:,end),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
% [X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,v(:,:,end),u(:,:,end),3,'black');
hold off;
subplot(1,2,2);
contourf(omz(:,:,end),'LevelStep',0.01,'EdgeColor','none');
colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;

%% wavenumber figure;

npad = 1000; % # of trailing zeros psi is to be padded with
             % for fft calculation
f = 8*((-npad/2):(npad/2))/(npad+1);
f = f(1:end-1);
f = fftshift(f);
[Xpsi,Ypsi] = meshgrid(f,f);
runsumYset = zeros(size(fft2(psi(:,:,end),npad,npad)));
starttimestep = trtu - round(3*trtu/4);
endtimestep = trtu;
for i = starttimestep:endtimestep
    Yset = fft2(psi(:,:,i),npad,npad);
    runsumYset = runsumYset + Yset;
end
Yse = runsumYset./(endtimestep-starttimestep);
Yse = abs(Yse/sqrt(npad*npad)).^2;
%%
figure;
contourf(Xpsi,Ypsi,Yse,'LevelStep',0.01,'EdgeColor','none');
colorbar; colormap jet;
axis square;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
box("on");
xlabel('$k_x$','Interpreter','latex','FontSize',30);
ylabel('$k_y$','Interpreter','latex','FontSize',30);

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
