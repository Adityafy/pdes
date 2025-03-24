%% Script for GSH Tangent Space with Spectral
% Using ETD1 explicit (predictor corrector)

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
% dt_ts = 0.1;
% gamma = 15;
trtu = 2000; % transient time units
% tstu = 10; % tangent space time units
seed = 3;
% tN = 1; % renormalization time units
% nnorm = tN/dt_ts;
% nnorm = 5;


Lx = 40;
Ly = Lx;
Nx = Lx;
Ny = Nx;
N = Nx*Ny;

% totmax = totu/dt_fpv;
nmax_tr = trtu/dt_tr;
% nmax_ts = tstu/dt_ts;

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
% kx(Nx/2+1) = sqrt(eps);
% ky(Ny/2+1) = sqrt(eps);
kx = fftshift(kx);
ky = fftshift(ky);
[Kx,Ky] = meshgrid(kx,ky);

L1 = epsilon - 1 - Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 ;
L2 = - sig * (Kx.^2 + Ky.^2) - sig * c^2;

%% Intervals
TrIntervals = floor(trtu);
% TsIntervals = floor(tstu);
% if tstu < 1
%     TsIntervals = 10;
% end


run_name = join(['specTr' num2str(Nx) 'eps' pointToDashed(epsilon) ...
            'sig' pointToDashed(sig) 'csq' pointToDashed(csq) ...
            'gm' pointToDashed(gm) 'trt' pointToDashed(trtu) ...
            's' pointToDashed(seed)]);

fprintf(join(['Running transients ' run_name ' ...\n']));
%% INITIAL CONDITIONS

% -----periodic ICs-----
% psitr(:,:,1) = 0.1*cos(X);
% psitrmat = psitr(:,:,1);
% % changing to vector
% psitrvec = latToVec(psitr);
% % assigning psihat for fft of psi
% psitrvec = fftshift(psitrvec);
% psih = fft(psitrvec); % psihat
% 
% zetatrmat = zeros(length(X), length(Y),1); % scalar field
% % zetatrmat = cos(X/16).*(1+sin(Y/16));
% zetatr = zeros(size(zetatrmat));
% utrmat = zeros(length(X), length(Y),1);
% vtrmat = zeros(length(X), length(Y),1);
% utr = zeros(length(X), length(Y),1);
% vtr = zeros(length(X), length(Y),1);
% zetatrvec = latToVec(zetatrmat);
% zetatrvec = fftshift(zetatrvec);
% zetatrh = fft(zetatrvec);
% 
% omztr = zeros(size(X));
% omztrmat = omztr;
% omztrvec = latToVec(omztr);
% omztrvec = fftshift(omztrvec);
% omztrh = fft(omztrvec);

% -----random ICs-----
rng(seed,"twister");
psitr(:,:,1) = -0.1 + (0.1+0.1) * rand(length(X), length(Y));
psitrmat = -0.1 + (0.1+0.1) * rand(length(X), length(Y));
% changing to vector
psitrvec = latToVec(psitr);
% assigning psihat for fft of psi
psitrvec = fftshift(psitrvec);
psih = fft(psitrvec); % psihat

%-------------- random IC for zeta and omega -----------------
zetatrmat = zeros(length(X), length(Y),1); % scalar field
% zetatrmat = cos(X/16).*(1+sin(Y/16));
zetatr = zeros(size(zetatrmat));
utrmat = zeros(length(X), length(Y),1);
vtrmat = zeros(length(X), length(Y),1);
utr = zeros(length(X), length(Y),1);
vtr = zeros(length(X), length(Y),1);
zetatrvec = latToVec(zetatrmat);
zetatrvec = fftshift(zetatrvec);
zetatrh = fft(zetatrvec);

rng(seed + 1,"twister");
omztr = -0.1 + (0.1+0.1)*rand(length(X), length(Y));
omztrmat = omztr;
omztrvec = latToVec(omztr);
omztrvec = fftshift(omztrvec);
omztrh = fft(omztrvec);

%% ETD explicit
Kdiff = Kx.^2 + Ky.^2;
firstrowKdiff = latToVec(Kdiff)';
KdiffMat = toeplitz(firstrowKdiff);
qx = latToVec(Kx);
qy = latToVec(Ky);
expL1hvec = exp(latToVec(L1*dt_tr));
expL2hvec = exp(latToVec(L2*dt_tr));
% KdiffMat = toeplitz(qx.^2+qy.^2);

tic;

for n = 1:nmax_tr
   
    % %------------------ zeta,  pseudospectral ------------------------
    zetamat2h = fft2(omztrmat)./(Kdiff+eps);
    zetatrh = fft(real(latToVec(ifft2(zetamat2h))));

    psihguess = psih;
    omzhguess = omztrh;
    for iters = 2
        psihpred = fft(real(ifft(psih))) .* expL1hvec ...
                    + 0.5 * dt_tr * N1hat(psihguess,zetatrh,qx,qy) ...
                    + 0.5 * dt_tr * expL1hvec .* N1hat(psih,zetatrh,qx,qy);
        omzhpred = fft(real(ifft(omztrh))).*expL2hvec ...
                    + 0.5 * dt_tr * N2hat(gm,psihguess,qx,qy) ...
                    + 0.5 * dt_tr * expL2hvec.* N2hat(gm,psih,qx,qy);
        psihguess = psihpred;
    end
    % -------------corrector, psi -------------
    psih = fft(real(ifft(psih))) .* expL1hvec ...
                    + 0.5 * dt_tr * N1hat(psihpred,zetatrh,qx,qy) ...
                    + 0.5 * dt_tr * expL1hvec .* N1hat(psih,zetatrh,qx,qy);

    % -------------corrector, omz -------------
    omztrh = fft(real(ifft(omztrh))) .* expL2hvec ...
                    + 0.5 * dt_tr * N2hat(gm,psihpred,qx,qy) ...
                    + 0.5 * dt_tr * expL2hvec.* N2hat(gm,psih,qx,qy);
    

    psitrmat = vecToLat(real(ifft(psih)),Nx,Ny);
    omztrmat = vecToLat(real(ifft(omztrh)),Nx,Ny);
    zetatrmat = vecToLat(real(ifft(zetatrh)),Nx,Ny);
    
    utrmat = vecToLat(real(ifft(1i*qy.*zetatrh)),Nx,Ny);
    vtrmat = vecToLat(real(ifft(-1i*qx.*zetatrh)),Nx,Ny);

    %------------- saving dynamics at intervals --------------
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

end

%%
figure; hold on;
attime = trtu;
contourf(X,Y,psitr(:,:,attime),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
% [X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,vtr(:,:,attime),utr(:,:,attime),3,'black');
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
quiver(X,Y,vtr(:,:,end),utr(:,:,end),3,'black');
hold off;
subplot(1,2,2);
contourf(omztr(:,:,end),'LevelStep',0.01,'EdgeColor','none');
colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;

%% saving
save(join([sfa run_name]), '-v7.3');

%% Nonlinear parts for GSH, transients and ts dynamics
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