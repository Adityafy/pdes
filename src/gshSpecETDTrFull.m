function [psitr,omztr,zetatr,utr,vtr] = gshSpecETDTrFull(cp,sv,trp)

epsilon = cp.epsilon;
sigma = cp.sig;
csq = cp.csq;
gm = cp.gm;

Lx = sv.Lx;
Ly = Lx;
Nx = Lx;
Ny = Nx;

dt = trp.dt;
trtu = trp.trtu;
seed = trp.seed;

nmax_tr = trtu/dt;

X = sv.X;
Y = sv.Y;
L1 = sv.L1;
L2 = sv.L2;

%% Intervals
TrIntervals = floor(trtu);

run_name = join(['specTr' num2str(Nx) 'eps' pointToDashed(epsilon) ...
            'sig' pointToDashed(sigma) 'csq' pointToDashed(csq) ...
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
Kdiff = sv.Kdiff;
qx = sv.qx;
qy = sv.qy;
expL1hvec = exp(latToVec(L1*dt));
expL2hvec = exp(latToVec(L2*dt));

tic;

for n = 1:nmax_tr
   
    % %------------------ zeta,  pseudospectral ------------------------
    zetamat2h = fft2(omztrmat)./(Kdiff+eps);
    zetatrh = fft(real(latToVec(ifft2(zetamat2h))));

    psihguess = psih;
    omzhguess = omztrh;
    for iters = 2
        psihpred = fft(real(ifft(psih))) .* expL1hvec ...
                    + 0.5 * dt * N1hat(psihguess,zetatrh,qx,qy) ...
                    + 0.5 * dt * expL1hvec .* N1hat(psih,zetatrh,qx,qy);
        omzhpred = fft(real(ifft(omztrh))).*expL2hvec ...
                    + 0.5 * dt * N2hat(gm,psihguess,qx,qy) ...
                    + 0.5 * dt * expL2hvec.* N2hat(gm,psih,qx,qy);
        psihguess = psihpred;
    end
    % -------------corrector, psi -------------
    psih = fft(real(ifft(psih))) .* expL1hvec ...
                    + 0.5 * dt * N1hat(psihpred,zetatrh,qx,qy) ...
                    + 0.5 * dt * expL1hvec .* N1hat(psih,zetatrh,qx,qy);

    % -------------corrector, omz -------------
    omztrh = fft(real(ifft(omztrh))) .* expL2hvec ...
                    + 0.5 * dt * N2hat(gm,psihpred,qx,qy) ...
                    + 0.5 * dt * expL2hvec.* N2hat(gm,psih,qx,qy);
    

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

end