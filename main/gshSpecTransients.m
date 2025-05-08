%% Script for GSH Tangent Space with Spectral
% Using ETD1 explicit (predictor corrector)
% With fft2 and ifft2
% No iterative used but functions are there if needed

clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% Parameter struct

% key parameters
epsilon = 0.7;
sig = 1;
csq = 0.1;
c = sqrt(csq);
gm = 50;
lam_0 = 2*pi; % critical roll wavelength
dt_tr = 0.1;
trtu = 2000; % transient time units
seed = 4;


nmax_tr = trtu/dt_tr;

wlmult = 8; % wavelength multiplier
Lx = lam_0 * wlmult;
Ly = Lx;
Nx = 8 * wlmult;
Ny = Nx;
N = Nx*Ny;

gamma = pi/2;
frgamma = gamma; %filtering radius gamma

% ---------  spatial grid in real domain -----------
x2 = linspace(-Lx/2,Lx/2,Nx+1)';
x = x2(1:Nx);
y2 = linspace(-Ly/2,Ly/2,Ny+1)';
y = y2(1:Ny);
dx = x(2) - x(1);
dy = y(2) - y(1);

[X,Y] = meshgrid(x,y);

%% Linear derivative operator for Fourier domain
% for both the equations.

kx = (2*pi/Lx)*(-Nx/2:Nx/2-1)';
ky = (2*pi/Ly)*(-Ny/2:Ny/2-1)';
kx = fftshift(kx);
ky = fftshift(ky);
% kx(1) = 10^-2;
% ky(1) = 10^-2;
% eps = 0;
[Kx,Ky] = meshgrid(kx,ky);

L1 = epsilon - 1 - Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 ;
L2 = - sig * (Kx.^2 + Ky.^2) - sig * c^2;

%% Intervals
TrIntervals = floor(trtu);
% TsIntervals = floor(tstu);
% if tstu < 1
%     TsIntervals = 10;
% end

run_name = join(['psetdL' num2str(ceil(Lx)) 'N' num2str(Nx) 'eps' pointToDashed(epsilon) ...
            'sig' pointToDashed(sig) 'csq' pointToDashed(csq) ...
            'gm' pointToDashed(gm) 'trt' pointToDashed(trtu) ...
            'fr' pointToDashed(frgamma) 's' pointToDashed(seed)]);

fprintf(join(['Running transients ' run_name ' ...\n']));

%% INITIAL CONDITIONS
% make sure zeta (zetatrmat) initialization is not too high (max < 0.1) or
% the code breaks.

%%%======= INITIALIZATION of Matrices =======
% subsections should be  commented/uncommented for periodic or random ICs

% -----periodic ICs-----
%%% note: the angular frequency should be an integral multiple of 2*pi/L
psitrmat = 0.01*(cos(6*pi*X/Lx)+sin(6*pi*Y/Ly));
zetatrmat = 0.01*(cos(2*pi*X/Lx)+sin(2*pi*Y/Ly));

% -----random ICs-----
% rng(seed,"twister");
% psitrmat = -0.1 + (0.1+0.1) * rand(length(X), length(Y));
% rng(seed + 1,"twister");
% zetatrmat = -0.1 + (0.1+0.1)*rand(length(X), length(Y));


%%%======= Initializing variables for time marching ========
psitr(:,:,1) = psitrmat; % 3D mat that's going to be appended and saved
psitrvec = latToVec(psitr); % changing to vector
% psitrvec = fftshift(psitrvec);
psihmat = fft2(psitrmat); % assigning psihat for fft of psi

zetatr = zetatrmat;
utrmat = zeros(length(X), length(Y),1);
vtrmat = zeros(length(X), length(Y),1);
utr = zeros(length(X), length(Y),1);
vtr = zeros(length(X), length(Y),1);
zetatrvec = latToVec(zetatrmat);
% zetatrvec = fftshift(zetatrvec);
zetatrhmat = fft2(zetatr);

% omztr = -fdlaplacian(zetatr, Nx,Ny,dx); % if fd laplacian is to be used
omztrmat = (Kx.^2+Ky.^2).*(zetatrmat);
omztr = omztrmat;
omztrhmat = fft2(omztrmat);

%% For FD gmres (if used)
%%% making the finite difference matrix A (fdA) for Ax=b
delsq = dx^2;
fdA = zeros(N);
I = 1:Nx;
J = 1:Ny;
Ip1 = circshift(I,-1);
Im1 = circshift(I,1);
Jp1 = circshift(J,-1);
Jm1 = circshift(J,1);
spotter = zeros(Nx,Ny);
iters = 1;
for i = 1:Nx
    for j = 1:Ny
        spotter(i,j) = 4/delsq;
        spotter(Ip1(i),J(j)) = -1/delsq;
        spotter(Im1(i),J(j)) = -1/delsq;
        spotter(I(i),Jp1(j)) = -1/delsq;
        spotter(I(i),Jm1(j)) = -1/delsq;
        fdA(iters,:) = latToVec(spotter)';
        iters = iters+1;
        spotter = zeros(Nx,Ny);
    end
end
fdA = sparse(fdA);

restart = 20;
tol = 1e-6; % error tolerance
maxit = 200; % max # of iterations

%% ETD explicit
Kdiff = Kx.^2 + Ky.^2;
firstrowKdiff = latToVec(Kdiff)';
KdiffMat = toeplitz(firstrowKdiff);
qx = latToVec(Kx);
qy = latToVec(Ky);
expL1hvec = exp(latToVec(L1*dt_tr));
expL2hvec = exp(latToVec(L2*dt_tr));

expL1hmat = vecToLat(expL1hvec,Nx,Ny);
expL2hmat = vecToLat(expL2hvec,Nx,Ny);
% KdiffMat = toeplitz(qx.^2+qy.^2);


%% Filtering
% fc = 0.25; % filtering constant
% gamma_filter = fc*pi;
% gammasq = gamma_filter^2;
% F_gamma = exp(-gammasq*(qx.^2+qy.^2)./2); % filtering on
% F_gamma = ones(size(qx)); % filtering off

%%
tic;
% figure;
% imagesc(omztrmat); colorbar; axis square; title('\Omega_z'); colormap jet;
    % subplot(1,2,1); imagesc(psitrmat); colorbar; axis square; title('\psi');
    % subplot(1,2,2); imagesc(omztrmat); colorbar; axis square; title('\Omega_z');
%     subplot(1,3,3); imagesc(zetatrmat); colorbar; axis square; title('\zeta');
Kdiff(1,1) = 1;
for n = 1:nmax_tr
    % close all; 
    % %================= zeta,  pseudospectral ======================
    
    % %----------------- spectral zeta calculation ------------------
    zetamat2h = fft2(omztrmat)./(Kdiff);
    zetamat2h(1,1) = 0+1i*0;
    zetamat2h(1,Ny) = 0+1i*0;
    zetamat2h(Nx,1) = 0+1i*0;
    zetamat2h(Nx,Ny) = 0+1i*0;
    zetatrhmat = fft2(real(ifft2(zetamat2h)));
    
    %----------------- gmres ------------------
    % zetamat2h = fft2(omztrmat)./(Kdiff+eps);
    % zetaguess = latToVec(real(ifft2(zetamat2h)));
    % zetaguess(1) = 0;
    % zetatrvec = gmres(fdA,latToVec(omztrmat),restart,tol,maxit, ...
    %                     [],[],zetatrvec);
    % zetatrh = fft(zetatrvec);

    %----------------- my iterative solver ------------------
    % zetatrvec = latToVec(iterZetaOmz(omztrmat,zetatrmat,dx,Nx,Ny));
    % zetatrh = fft(zetatrvec);
    
    %==================== psi, omz,  pseudospectral, etd ==================
    psihguess = psihmat;
    omzhguess = omztrhmat;
    for iters = 2
        psihpred = fft2(real(ifft2(psihmat))) .* expL1hmat ...
                    + 0.5 * dt_tr * N1hat(psihguess,zetatrhmat,Kx,Ky) ...
                    + 0.5 * dt_tr * expL1hmat .* N1hat(psihmat,zetatrhmat,Kx,Ky);
        % omzhpred = fft(real(ifft(omztrh))).*expL2hvec ...
                    % + 0.5 * dt_tr * F_gamma.*N2hat(gm,psihguess,qx,qy) ...
                    % + 0.5 * dt_tr * expL2hvec.* F_gamma.*N2hat(gm,psih,qx,qy);
        omzhpred = fft2(real(ifft2(omztrhmat))).*expL2hmat ...
                    + 0.5 * dt_tr * N2hat(gamma,gm,psihguess,Kx,Ky) ...
                    + 0.5 * dt_tr * expL2hmat.*N2hat(gamma,gm,psihmat,Kx,Ky);
        psihguess = psihpred;
        omzhguess = omzhpred;
    end
    % -------------corrector, psi -------------
    psihmat = fft2(real(ifft2(psihmat))) .* expL1hmat ...
                    + 0.5 * dt_tr * N1hat(psihpred,zetatrhmat,Kx,Ky) ...
                    + 0.5 * dt_tr * expL1hmat .* N1hat(psihmat,zetatrhmat,Kx,Ky);

    % -------------corrector, omz -------------
    % omztrh = fft(real(ifft(omztrh))) .* expL2hvec ...
    %                 + 0.5 * dt_tr * F_gamma.*N2hat(gm,psihpred,qx,qy) ...
    %                 + 0.5 * dt_tr * expL2hvec.* F_gamma.*N2hat(gm,psih,qx,qy);
    omztrhmat = fft2(real(ifft2(omztrhmat))) .* expL2hmat ...
                    + 0.5 * dt_tr * N2hat(gamma,gm,psihpred,Kx,Ky) ...
                    + 0.5 * dt_tr * expL2hmat.*N2hat(gamma,gm,psihmat,Kx,Ky);
    

    psitrmat = real(ifft2(psihmat,'symmetric'));
    omztrmat = real(ifft2(omztrhmat,'symmetric'));
    zetatrmat = real(ifft2(zetatrhmat,'symmetric'));
    
    utrmat = real(ifft2(1i*Ky.*zetatrhmat,'symmetric'));
    vtrmat = real(ifft2(-1i*Kx.*zetatrhmat,'symmetric'));

    %================== saving dynamics at intervals ==================
    if rem(n,nmax_tr/TrIntervals) == 0
        psitr(:,:,round(n*TrIntervals/nmax_tr)) = psitrmat;
        omztr(:,:,round(n*TrIntervals/nmax_tr)) = omztrmat;
        zetatr(:,:,round(n*TrIntervals/nmax_tr)) = zetatrmat;
        utr(:,:,round(n*TrIntervals/nmax_tr)) = utrmat;
        vtr(:,:,round(n*TrIntervals/nmax_tr)) = vtrmat;
        
    end
    if rem(n,round(nmax_tr/10)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax_tr);
        toc;
    end
    if isnan(abs(sum(sum(psitrmat))))
        error('blow up occured in psi..');
    end
    
    % subplot(1,2,1); 
    % imagesc([psitrmat psitrmat]); colorbar; colormap jet; clim([-0.8 0.8]);
    % subplot(1,2,2); 
    imagesc(psitrmat); colorbar; axis square; clim([-1 1]); colormap jet;
    % % subplot(1,3,3); imagesc(zetatrmat); colorbar; axis square;
    drawnow;
end

%%
figure; hold on;
attime = trtu;
contourf(X,Y,psitr(:,:,attime),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
% [X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,utr(:,:,attime),vtr(:,:,attime),3,'black');
box on;
xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
hold off;

%%
figure;
subplot(1,2,1);
hold on;
attime = trtu;
contourf(X,Y,psitr(:,:,end),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
% [X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,utr(:,:,end),vtr(:,:,end),3,'black');
hold off;
subplot(1,2,2);
contourf(omztr(:,:,end),'LevelStep',0.01,'EdgeColor','none');
colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;

%%
figure;
subplot(1,2,1);
hold on;
attime = trtu;
contourf(X,Y,psitr(:,:,end),'LevelStep',0.01,'EdgeColor','none');
xlim([-Lx/2 Lx/2-1]); ylim([-Ly/2 Ly/2-1]); title '\psi';
clim([-1 1]);
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
quiver(X,Y,utr(:,:,end),vtr(:,:,end),3,'black');
hold off;
subplot(1,2,2);
imagesc(x,y,omztr(:,:,end)); clim([-2.5 2.5]); title '\Omega_z';
colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;

%%
figure;
imagesc(x,y,omztr(:,:,end-3)); colorbar; colormap jet; axis square; clim([-2.5 2.5])
title '\Omega_z';

%% wavenumber calculation

npad = 1000; % # of trailing zeros psi is to be padded with
             % for fft calculation
knpad = wlmult*((-npad/2):(npad/2))/(npad+1);
knpad = knpad(1:end-1);
knpad = fftshift(knpad);
[Kxn,Kyn] = meshgrid(knpad,knpad);
runsumPsiHat = zeros(size(fft2(psitr(:,:,end),npad,npad)));
starttimestep = trtu - round(3*trtu/4);
endtimestep = trtu;
for time_step = starttimestep:endtimestep
    psihat_t = fft2(psitr(:,:,time_step),npad,npad);
    runsumPsiHat = runsumPsiHat + psihat_t;
end
psihat = runsumPsiHat./(endtimestep-starttimestep);
psihat = abs(psihat/sqrt(npad*npad)).^2;

% wavenumber figure
figure;
contourf(Kxn,Kyn,psihat,'LevelStep',0.0005,'EdgeColor','none');
colorbar; colormap jet;
axis square;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
box("on");
xlabel('$k_x$','Interpreter','latex','FontSize',30);
ylabel('$k_y$','Interpreter','latex','FontSize',30);

%% saving
save(join([sfa run_name]), '-v7.3');

%% Nonlinear parts for GSH, transients and ts dynamics
function N1h = N1hat(psih,zetah,Kx,Ky)
% function for the nonlinear part of equation 1
    N1h = -fft2(real((ifft2(psih)).^3)) ...
            - fft2(real( ifft2(1i*Ky.*zetah) .* ifft2(1i*Kx.*psih) )) ...
                + fft2(real( ifft2(1i*Kx.*zetah) .* ifft2(1i*Ky.*psih) ));
end

function N2h = N2hat(gamma,gm,psih,Kx,Ky)
% function for the nonlinear part of equation 2

    psi_y = ifft2(1i*Ky.*psih);
    psi_xxx = ifft2(-1i*Kx.^3.*psih);
    psi_xyy = ifft2(-1i*Kx.*Ky.^2.*psih);
    psi_x = ifft2(1i*Kx.*psih);
    psi_yxx = ifft2(-1i*Ky.*Kx.^2.*psih);
    psi_yyy = ifft2(-1i*Ky.^3.*psih);

    if gamma == 0
        % fft of the nonlinear part
        N2h =  fft2(real(-gm *  psi_y.*(psi_xxx+psi_xyy))) ...
            + fft2(real(gm * psi_x.*(psi_yxx+psi_yyy))) ;

    elseif gamma > 0
        %%% taking gaussian filter of rhs (seems to work)
        rhs_wo_gm = -psi_y.*(psi_xxx+psi_xyy) + psi_x.*(psi_yxx+psi_yyy);
        rhs_wo_gm_matrix = real(rhs_wo_gm);
        filtered_rhs_wo_gm = imgaussfilt(rhs_wo_gm_matrix,gamma, ...
            'FilterDomain','spatial','Padding','circular');
        filtered_rhs_matrix = gm * filtered_rhs_wo_gm;
        N2h = fft2(real(filtered_rhs_matrix));
        % imagesc(gm*[filtered_rhs_wo_gm filtered_rhs_wo_gm]); colorbar; colormap jet;
        % drawnow;
    end
    

    %%% fft of the nonlinear part if F_gamma is applied to the N2h (this
    %%% seems to be working in a weird way)
    % N2h =  F_gamma.*( fft(real(-gm *  psi_y.*(psi_xxx+psi_xyy))) ...
    %      + fft(real(gm * psi_x.*(psi_yxx+psi_yyy))) ) ;

    %%% if F_gamma is added multiplied within the fft brackets (doesn't
    %%% work)
    % N2h =  fft( F_gamma.*real(-gm *  psi_y.*(psi_xxx+psi_xyy)) ...
    %      + F_gamma.*real(gm * psi_x.*(psi_yxx+psi_yyy))) ;

    %%% if F_gamma is added to the derivatives (doesn't work)
    % psi_y = ifft(1i*F_gamma.*qy.*psih);
    % psi_xxx = ifft((-1i)*F_gamma.*qx.^3.*psih);
    % psi_xyy = ifft((-1i)*F_gamma.*qx.*qy.^2.*psih);
    % psi_x = ifft(1i*F_gamma.*qx.*psih);
    % psi_yxx = ifft((-1i)*F_gamma.*qy.*qx.^2.*psih);
    % psi_yyy = ifft((-1i)*F_gamma.*qy.^3.*psih);
end

function g = fdlaplacian(f, Nx, Ny, dx)
% zeta = fdlaplacian(omzmat, zetamat)
% Laplacian calculator using finite difference
% Computes g = \nabla^2 (f)
% with periodic boundary conditions.
    I = 1:Nx;
    J = 1:Ny;
    Ip1 = circshift(I,-1);
    Im1 = circshift(I,1);
    Jp1 = circshift(J,-1);
    Jm1 = circshift(J,1);
    g = zeros(size(f));
    for i = 1:length(I)
        for j = 1:length(J)
            g(i,j) = (1/(dx^2)) * ...
                ( f(Im1(i),J(j)) + f(Ip1(i),J(j)) ...
                - 4 * f(I(i),J(j)) ...
                + f(I(i),Jm1(j)) + f(I(i),Jp1(j)) );
        end
    end
end

function zeta = iterZetaOmz(omzmat, zetamat, dx,Nx,Ny)
% zeta = iterativeZetaOmz(omzmat, zetamat, p)
% An iterative poisson solver for
% Omega_z = - Laplacian(Zeta)
% with periodic boundary conditions.
% mind the negative sign if used for other kinds of Poisson solvers
% Number of iterations is fixed for now (= 800)
    delta = dx;

    I = 1:Nx;
    J = 1:Ny;
    Ip1 = circshift(I,-1);
    Im1 = circshift(I,1);
    Jp1 = circshift(J,-1);
    Jm1 = circshift(J,1);

    zeta = zetamat;
    iterations = 500;
    difference = zeros(1,iterations);

    for k = 1:iterations
        for i = 1:length(I)
            for j = 1:length(J)
                zeta(i,j) = ((delta^2)/4) * omzmat(I(i),J(j)) ...
                    + (1/4) * ( zetamat(Im1(i),J(j)) + zetamat(Ip1(i),J(j)) ...
                    + zetamat(I(i),Jm1(j)) + zetamat(I(i),Jp1(j)) );
            end
        end

        lapl = zeros(size(zeta));
        for i = 1:length(I)
            for j = 1:length(J)
                lapl(i,j) = (1/(delta^2)) * ...
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