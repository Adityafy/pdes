%% Script for GSH Tangent Space with Spectral
% Using ETD1 explicit (predictor corrector)
% transient file must be loaded.

clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% %%%%%%%% Parameters
dt_ts = 0.1;
tstu = 100; % tangent space time units
tN = 1; % renormalization time units
nnorm = tN/dt_ts;
nmax_ts = tstu/dt_ts;
% nmax_ts = 50000;
TsIntervals = floor(tstu);
% TsIntervals = 100;
if tstu < 1
    TsIntervals = 10;
end

%% Load transients
load('../../pdesDataDump/specTr64eps0-7sig1csq4gm0trt1000s3.mat');

%% Dynamics after transients initialization

% -----random ICs-----
rng(seed,"twister");
psi(:,:,1) = psitr(:,:,end);
psimat = psi(:,:,1);
psivec = latToVec(psi); % changing to vector
psivec = fftshift(psivec); % assigning psihat for fft of psi
psih = fft(psivec); % psihat

%-------------- random IC for zeta and omega -----------------
zetamat = zetatrmat; % scalar field
zeta = zetatr(:,:,end);
umat = utrmat;
vmat = vtrmat;
u = utr(:,:,end);
v = vtr(:,:,end);
zetavec = latToVec(zetamat);
zetavec = fftshift(zetavec);
zetah = fft(zetavec);

omz = omztr(:,:,end);
omzmat = omz;
omzvec = latToVec(omz);
omzvec = fftshift(omzvec);
omzh = fft(omzvec);

%% Tangent space initialization
nv = 1; % number of vectors to be calculated nv <= N
rng(seed + 2, "twister");
% dpsi = -0.1 + (0.1+0.1) * rand(N,nv);
% 
% rng(seed + 3, "twister");
% domz = -0.1 + (0.1+0.1) * rand(N,nv);

dpsih = zeros(N,nv);
domzh = zeros(N,nv);
dH = rand(2*N,nv,1); % columns of dH are perturbation vectors
dH = 0.001*orth(dH); % make them orthonormal and small
dpsi = dH(1:N,:);
domz = dH(N+1:2*N,:);
for iter = 1:nv
    dpsih(:,iter) = fft(dpsi(:,iter));
    domzh(:,iter) = fft(domz(:,iter));
end
dH(:,:,1) = [dpsi; domz]; % initializing perturbation vectors

laminst = zeros(nv,1); % lambda(k) instantaneous

dHmag = zeros(nv,nmax_ts); % magnitude of the perturbation vectors

dH1 = zeros(2*N,nmax_ts/nnorm);

%% ETD explicit
Kdiff = Kx.^2 + Ky.^2;
firstrowKdiff = latToVec(Kdiff)';
KdiffMat = toeplitz(firstrowKdiff);
qx = latToVec(Kx);
qy = latToVec(Ky);
expL1hvec = exp(latToVec(L1*dt_ts));
expL2hvec = exp(latToVec(L2*dt_ts));
% KdiffMat = toeplitz(qx.^2+qy.^2);
L1vec = epsilon - 1 - qx.^4 - 2*(qx.^2).*(qy.^2) ...
            - qy.^4 + 2*qx.^2 + 2*qy.^2 ;
tic;

for n = 1:nmax_ts
   
    % %------------------ zeta,  pseudospectral ------------------------
    zetamat2h = fft2(omzmat)./(Kdiff+eps);
    zetah = fft(real(latToVec(ifft2(zetamat2h))));
    % zetamat = ifft2(zetamat2h);
    % poisson_test = isequal(omzmat,ifft2(Kdiff*zetamat2h));
    % zetah = fft(real(latToVec(ifft2(zeta2h))));

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
                    + 0.5 * dt_ts * N1hat(psihguess,zetah,qx,qy) ...
                    + 0.5 * dt_ts * expL1hvec .* N1hat(psih,zetah,qx,qy);
        omzhpred = fft(real(ifft(omzh))).*expL2hvec ...
                    + 0.5 * dt_ts * N2hat(gm,psihguess,qx,qy) ...
                    + 0.5 * dt_ts * expL2hvec.* N2hat(gm,psih,qx,qy);
        psihguess = psihpred;
    end
    % -------------corrector, psi -------------
    psih = fft(real(ifft(psih))) .* expL1hvec ...
                    + 0.5 * dt_ts * N1hat(psihpred,zetah,qx,qy) ...
                    + 0.5 * dt_ts * expL1hvec .* N1hat(psih,zetah,qx,qy);

    % -------------corrector, omz -------------
    omzh = fft(real(ifft(omzh))) .* expL2hvec ...
                    + 0.5 * dt_ts * N2hat(gm,psihpred,qx,qy) ...
                    + 0.5 * dt_ts * expL2hvec.* N2hat(gm,psih,qx,qy);
    

    psimat = vecToLat(real(ifft(psih)),Nx,Ny);
    omzmat = vecToLat(real(ifft(omzh)),Nx,Ny);
    zetamat = vecToLat(real(ifft(zetah)),Nx,Ny);
    
    umat = vecToLat(real(ifft(1i*qy.*zetah)),Nx,Ny);
    vmat = vecToLat(real(ifft(-1i*qx.*zetah)),Nx,Ny);

    %===================== tangent space ========================
    %--------------------- coefficients ----------------------------
    alpha1 = coeffTS1dpsi(L1vec,psimat,psih,zetah,qx,qy);
    % alpha1 = coeffTS1dpsi(L1,psimat,omzh,qx,qy);
    beta1 = coeffTS1domz(psih,qx,qy);
    alpha2 = coeffTS2dpsi(gm,psih,qx,qy);
    beta2 = coeffTS2domz(sig,csq,qx,qy);
    %--------------------- ETD (method1) ----------------------------
    dpsih_temp = dpsih;
    domzh_temp = domzh;
    for k = 1:nv
        dpsih(:,k) = fft(real(ifft(dpsih_temp(:,k)))) .* exp(alpha1*dt_ts) ...
            + fft(real(ifft(domzh_temp(:,k)))) .* exp(beta1*dt_ts);
        domzh(:,k) = fft(real(ifft(dpsih_temp(:,k)))) .* exp(alpha2*dt_ts) ...
            + fft(real(ifft(domzh_temp(:,k)))) .* exp(beta2*dt_ts);
        dH(:,k) = [real(ifft(dpsih(:,k))); real(ifft(domzh(:,k)))];
    end
    %--------------------- ETD (method2) ----------------------------
    % alpha1 = c1dpsi(L1,psimat,zetamat,Kx,Ky);
    % beta1 = c1domz(psimat,Kx,Ky);
    % alpha2 = c2dpsi(gm,psimat,Kx,Ky);
    % beta2 = c2domz(L2);
    % expalpha1h = exp(alpha1*dt_ts);
    % expbeta1h = exp(beta1*dt_ts);
    % expalpha2h = exp(alpha2*dt_ts);
    % expbeta2h = exp(beta2*dt_ts);
    % for k = 1:nv
    %     dpsih(:,k) = fft(real(ifft(dpsih(:,k)))) .* expalpha1h ...
    %         + fft(real(ifft(domzh(:,k)))) .* expbeta1h;
    %     domzh(:,k) = fft(real(ifft(dpsih(:,k)))) .* expalpha2h ...
    %         + fft(real(ifft(domzh(:,k)))) .* expbeta2h;
    %     dH(:,k) = [real(ifft(dpsih(:,k))); real(ifft(domzh(:,k)))];
    % end
    %--------------------------- FOE ----------------------------
    % for k = 1:nv
    %     dpsih(:,k) = (alpha1.*dt_ts + 1) .* dpsih(:,k) ...
    %                     + beta1 .* dt_ts .* domzh(:,k);
    %     domzh(:,k) = (alpha2.*dt_ts) .* dpsih(:,k) ...
    %                     + (beta2 .* dt_ts + 1) .* domzh(:,k);
    %     dH(:,k) = [real(ifft(dpsih(:,k))); real(ifft(domzh(:,k)))];
    % end

    % --------- renormalization and LLE ----------
    if rem(n,nnorm) == 0
        [Q,R]= qr(dH);
        dH = Q(:,1:nv);
        dH1(:,n/nnorm) = dH(:,1);
        % laminst(:,n/nnorm) = (1/tN) * log(abs(diag(R)));
        for k = 1:nv
            laminst(k,n/nnorm) = (1/tN)*log(abs(R(k,k)));
        end
    end
    % extracting dpsih and domzh from dH and calculating ||dH(k)||
    for k = 1:nv
        dpsih(:,k) = fft(real(dH(1:N,k)));
        domzh(:,k) = fft(real(dH(N+1:2*N,k)));
        dHmag(k,n) = norm(dH(:,k));
    end

    % --------- separate renormalization and LLE ----------
    % if rem(n,nnorm) == 0
    %     [Qdpsi,Rdpsi]= qr(dH(1:N,1:nv));
    %     dpsi_ortho = Qdpsi(:,1:nv);
    %     dpsi = dpsi_ortho;
    %     % dH1(:,n/nnorm) = dpsi_ortho(:,1);
    %     % laminst(:,n/nnorm) = (1/tN) * log(abs(diag(R)));
    %     for k = 1:nv
    %         laminst(k,n/nnorm) = (1/tN)*log(abs(Rdpsi(k,k)));
    %     end
    %     [Qomz,~]= qr(dH(N+1:2*N,1:nv));
    %     domz_ortho = Qomz(:,1:nv);
    %     domz = domz_ortho;
    %     % extracting dpsih and domzh from dpsi and domz
    %     dH = [dpsi; domz];
    %     for k = 1:nv
    %         dpsih(:,k) = fft(dpsi(:,k));
    %         domzh(:,k) = fft(domz(:,k));
    %         dH(:,k) = dH(:,k)./norm(dH(:,k));
    %         dHmag(k,n) = norm(dH(:,k));
    %     end
    % 
    % end


    %------------- saving dynamics at intervals --------------
    if rem(n,nmax_ts/TsIntervals) == 0
        psi(:,:,round(n*TsIntervals/nmax_ts)) = psimat;
        omz(:,:,round(n*TsIntervals/nmax_ts)) = omzmat;
        zeta(:,:,round(n*TsIntervals/nmax_ts)) = zetamat;
        u(:,:,round(n*TsIntervals/nmax_ts)) = umat;
        v(:,:,round(n*TsIntervals/nmax_ts)) = vmat;
        pertvecs(:,:,round(n*TsIntervals/nmax_ts)) = dH;
    end
    if rem(n,round(nmax_ts/10)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax_ts);
        toc;
    end
    if isnan(abs(sum(sum(psimat))))
        error('blow up occured in psi..');
    end
    if n>nmax_ts/3
        csq = 0.5;
        c = sqrt(csq);
        gm = 50;
    end

end

%%
lam1gs = (1/(tstu/tN))*mean(laminst);
%%
% figure; contourf(psi(:,:,4),'LevelStep',0.01,'EdgeColor','none');
% colormap jet; colorbar;
% set(gca,'YDir','normal'); axis square;

%%
figure; hold on;
attime = tstu;
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
plot1 = contourf(psi(:,:,tstu),'levels',0.1, 'Linecolor', 'none');
set(gca,'YDir','normal');
hold on;
colormap jet;
clim([-0.8 0.8]);
colorbar(gca,'Ticks',[-0.8 0.8]);
plot2 = quiver(X,Y,v(:,:,tstu),u(:,:,tstu),3,'black');
xlim([1 Nx]);
ylim([1 Ny]);
axis square;
%%
figure;
subplot(1,2,1);
hold on;
attime = tstu;
contourf(X,Y,psi(:,:,end),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
% [X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,v(:,:,end),u(:,:,end),3,'black');
xlim([X(1) X(end)]);
ylim([Y(1) Y(end)]);
hold off;
box on;
subtitle('$\psi$','Interpreter','latex');
subplot(1,2,2);
contourf(omz(:,:,end),'LevelStep',0.01,'EdgeColor','none');
colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;
subtitle('$\Omega z$','Interpreter','latex');

%% psi at end
figure;
hold on;
attime = tstu;
contourf(X,Y,psi(:,:,end),'LevelStep',0.01,'EdgeColor','none');
set(gca,'YDir','normal'); axis square;
colormap jet; colorbar;
% [X,Y] = meshgrid(1:Nx,1:Ny);
quiver(X,Y,v(:,:,end),u(:,:,end),3,'black');
xlim([X(1) X(end)]);
ylim([Y(1) Y(end)]);
hold off;
box on;
subtitle('$\psi$','Interpreter','latex', 'FontSize',20);

%% omz
figure;
contourf(omz(:,:,end),'LevelStep',0.01,'EdgeColor','none');
colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;
subtitle('$\Omega z$','Interpreter','latex');

%%
dpsi1 = vecToLat(real(pertvecs(1:N,1,end)),Nx,Ny);
figure; imagesc(dpsi1'); colormap jet; colorbar;
set(gca,'YDir','normal');
set(gca,'Box','on');
axis square;
subtitle('$\delta \psi_1$','Interpreter','latex');

%%
figure; hold on;
dpsi1 = vecToLat(real(pertvecs(1:N,1,end)),Nx,Ny);

contourf(abs(dpsi1'),'levels',0.0001, 'Linecolor', 'none'); colormap jet; colorbar;
contourf(psi(:,:,end),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 2);
set(gca,'YDir','normal');
set(gca,'Box','on');
axis square;
subtitle('$|\delta \psi_1|$','Interpreter','latex','FontSize',20);

%%
figure; plot(dHmag,'-o','LineWidth',1);
set(gca,'TickLabelInterpreter','tex','FontSize',15);
xlabel('$n$', 'FontSize',30,'Interpreter','latex');
ylabel('$|d\vec{H}^{(1)}|$', 'FontSize',30, 'Interpreter','latex','Rotation',0);
set(gca,'Box', 'on', 'LineWidth',1.5);
axis square;
%%
figure; loglog(dHmag, '-o','LineWidth',1);
set(gca,'TickLabelInterpreter','tex','FontSize',15);
xlabel('$n$', 'FontSize',30,'Interpreter','latex');
ylabel('$|d\vec{H}^{(1)}|$', 'FontSize',30, 'Interpreter','latex','Rotation',0);
set(gca,'Box', 'on', 'LineWidth',1.5);
axis square;
%%
figure; plot(laminst(1,:), '-o','LineWidth',1.5);
set(gca,'TickLabelInterpreter','tex','FontSize',15);
xlabel('$n$', 'FontSize',30,'Interpreter','latex');
ylabel('$\lambda_1(t)$', 'FontSize',30, 'Interpreter','latex','Rotation',0);
set(gca,'Box', 'on', 'LineWidth',1.5);
axis square;

%% lambda_1 cummulative average (lam1ca)
lam1cafig = figure;
axlam1fig = axes('Parent',lam1cafig);
hold(axlam1fig,'on');
lam1inst = laminst(1,:);

time = linspace(0,tstu,floor(tstu/tN));
lam1ca = cumsum(lam1inst)./(1:length(lam1inst));
plot(time,lam1ca,'MarkerSize',10,'Marker','.','LineWidth',1.5);
% plot(cumsum(lam1inst)./(1:length(lam1inst)), '-o');
xlim([0 time(end)]);
yline(0,'Parent',axlam1fig);

set(gca, 'box', 'on');
axis(axlam1fig,'square');
hold(axlam1fig,'off');
set(axlam1fig,'FontSize',15,'LineWidth',1.5,'XMinorTick','on','YMinorTick','on');
ylabel('$\left<\lambda_{1}(t)\right>$','FontSize',30,'Interpreter','latex','Rotation',0);
xlabel('$t$','FontSize',30,'Interpreter','latex');
%% 
for iter = 1:length(pertvecs(1,1,:))
    dpsi1mat(:,:,iter) = vecToLat(real(pertvecs(1:N,1,iter)),Nx,Ny);
end
%% wavenumber figure;

npad = 1000; % # of trailing zeros psi is to be padded with
             % for fft calculation
f = 8*((-npad/2):(npad/2))/(npad+1);
f = f(1:end-1);
f = fftshift(f);
[Xpsi,Ypsi] = meshgrid(f,f);
runsumYset = zeros(size(fft2(psi(:,:,end),npad,npad)));
starttimestep = tstu - round(3*tstu/4);
endtimestep = tstu;
for i = starttimestep:endtimestep
    Yset = fft2(psi(:,:,i),npad,npad);
    runsumYset = runsumYset + Yset;
end
Yse = runsumYset./(endtimestep-starttimestep);
Yse = abs(Yse/sqrt(npad*npad)).^2;

%%
figure;
contourf(Xpsi,Ypsi,Yse,'LevelStep',0.0005,'EdgeColor','none');
colorbar; colormap jet;
axis square;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
box("on");
xlabel('$k_x$','Interpreter','latex','FontSize',30);
ylabel('$k_y$','Interpreter','latex','FontSize',30);


%% functions

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

%% Functions for TS coefficients, vectorization first calculation later
function Lvec = coeffTS1dpsi(L1vec,psimat,psih,zetah,qx,qy)
% coefficient of \vec{\delta\psi}^{(k)} in the first TS equation
Lvec = L1vec - 3.*(real(ifft(psih).^2)) - (1i.*qx) .* (real(ifft(1i*qy.*zetah))) ...
    +  (1i.*qy) .* (real(ifft(1i*qx.*zetah))) ;
end

function Lvec = coeffTS1domz(psih,qx,qy)
% coefficient of \vec{\delta\Omega_z}^{(k)} in the first TS equation
Lvec = - real(ifft(1i*qx.*psih)) .* (1i*qy./(qx.^2+qy.^2+eps)) ...
         + real(ifft(1i*qy.*psih)) .* (1i*qx./(qx.^2+qy.^2+eps));

end

function Lvec = coeffTS2dpsi(gm,psih,qx,qy)
% coefficient of \vec{\delta\psi}^{(k)} in the second TS equation
Lvec = gm * real(ifft(1i*qy.*psih)) .* (1i*qx.^3 + 1i*qx.*qy.^2) ...
        -gm * real(ifft(1i*qx.*psih)) .* (1i*qx.^2.*qy + 1i*qy.^3) ...
        - gm * 1i*qy .* (real(ifft(-1i*qx.^3.*psih))+real(ifft(-1i*qx.*qy.^2.*psih))) ...
        + gm * 1i*qx .* (real(ifft(-1i*qy.*qx.^2.*psih))+real(ifft(-1i*qy.^3.*psih)));
end

function Lvec = coeffTS2domz(sig,csq,qx,qy)
% coefficient of \vec{\delta\Omega_z}^{(k)} in the second TS equation
Lvec = - sig*(qx.^2+qy.^2) - sig*csq;
end

%% Functions for TS coefficients -- calculation first, vectorization later
% function alpha1 = c1dpsi(L1,psimat,zetamat,Kx,Ky)
% alpha1 = L1 - 3*psimat.^2 - (1i.*Kx) .* real(ifft2(1i*Ky.*fft2(zetamat))) ...
%     +  (1i.*Ky) .* real(ifft2(1i*Kx.*fft2(zetamat))) ;
% alpha1 = latToVec(alpha1);
% % expalpha1h = exp(latToVec(alpha1*dt_ts+eye(size(alpha1))));
% end
% 
% function beta1 = c1domz(psimat,Kx,Ky)
% beta1 = - real(ifft2(1i*Kx.*fft2(psimat))) .* (1i*Ky./(Kx.^2+Ky.^2+eps)) ...
%          + real(ifft2(1i*Ky.*fft2(psimat))) .* (1i*Kx./(Kx.^2+Ky.^2+eps));
% beta1 = latToVec(beta1);
% % expbeta1h = exp(latToVec(beta1)*dt_ts);
% % firstrow = Lvec';
% % coeffmat2 = toeplitz(firstrow);
% end
% 
% function alpha2 = c2dpsi(gm,psimat,Kx,Ky)
% alpha2 = gm * real(ifft2(1i*Ky.*fft2(psimat))) .* (1i*Kx.^3 + 1i*Kx.*Ky.^2) ...
%         -gm * real(ifft2(1i*Kx.*fft2(psimat))) .* (1i*Kx.^2.*Ky + 1i*Ky.^3) ...
%         + gm * 1i*Ky .* (real(ifft2(1i*Kx.^3.*fft2(psimat)))+ ...
%                             real(ifft2(1i*Kx.*Ky.^2.*fft2(psimat)))) ...
%         - gm * 1i*Kx .* (real(ifft2(1i*Ky.*Kx.^2.*fft2(psimat))) ...
%                             +real(ifft2(1i*Ky.^3.*fft2(psimat))));
% alpha2 = latToVec(alpha2);
% % expalpha2h = exp(latToVec(alpha2)*dt_ts);
% end
% 
% function beta2 = c2domz(L2)
% beta2 = latToVec(L2);
% % beta2 = latToVec(L2*dt_ts+eye(size(L2)));
% % firstrow = Lvec';
% % coeffmat2 = toeplitz(firstrow);
% end 