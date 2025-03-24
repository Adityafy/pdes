clear all;
close all;
addpath('../src/');

mfa = '../../pdesDataDump/media/'; % media folder address
dfa = '../../pdesDataDump/data/'; % saved data folder address

%% Parameters

dt_tr = 0.1;
% dt_ts = 0.1;
dt_ts = 1e-6;
trtu = 2000; % transient time units
tstu = 2; % tangent space time units
seed = 4;
tN = 0.5; % renormalization time units
nnorm = tN/dt_ts;

Nx = 64;

nmax_tr = trtu/dt_tr;
nmax_ts = tstu/dt_ts;

% ---------  spatial grid in real domain -----------
x2 = (pi)*linspace(-Nx/2,Nx/2,Nx+1)';
x = x2(1:Nx);
% x = 32*pi*(1:Nx)'/Nx;
dx = x(2) - x(1);

% wavenumber grid
kx = (2*pi/Nx)*(-Nx/2:Nx/2-1)';
% kx = [0:Nx/2-1 0 -Nx/2+1:-1]'/16;
kx = fftshift(kx);

alpha = kx.^2 - kx.^4;

% -----periodic ICs-----
utr = cos(x/16).*(1+sin(x/16)); 
nn = 0;

% -----random ICs-----
% rng(seed,"twister");
% utr = zeros(size(x));
% utr(:,1) = 0.01*rand(size(utr));

% -----In Fourier space----
utrhat = fft(utr(:,1)); % initializing Fourier space variable

expalphadt = exp(alpha * dt_tr);

tic
for n = 1:nmax_tr
    % predictor-corrector ETD:
    uguess = utrhat;
    for iters = 4
        upred = fft(real(ifft(utrhat))) .* expalphadt ...
                    + 0.5 * dt_tr * N1hat(uguess,kx) ...
                    + 0.5 * dt_tr * expalphadt .* N1hat(utrhat,kx);
        uguess = upred;
    end
    utrhat = fft(real(ifft(utrhat))) .* expalphadt ...
                    + 0.5 * dt_tr * N1hat(upred,kx) ...
                    + 0.5 * dt_tr * expalphadt .* N1hat(utrhat,kx);
    
    utr(:,n+1) = real(ifft(utrhat));
    nn = [nn n];
    if rem(n,round(nmax_tr/10)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax_tr);
        toc;
    end
    if isnan(abs(sum(utr)))
        error('blow up occured in u ..');
    end

end
toc
%%
figure; imagesc(utr');
colormap jet; colorbar;
clim([-4 4]);
set(gca,"YDir","normal");

%%
% figure;
% surf(nn,x,utr), shading interp, lighting phong, axis tight 
% view([-90 90]), colormap(autumn); set(gca,'zlim',[-5 50]) 
% light('color',[1 1 0],'position',[-1,2,2]) 
% material([0.30 0.60 0.60 40.00 1.00]);

%% Tangent space

% dynamics:
uts = real(ifft(utrhat));
u(:,1) = uts; 
nnts = 0;
utshat = utrhat; % initializing Fourier space variable
expalphadtts = exp(alpha * dt_tr);

% tangent space variables:
nv = 1; % number of vectors to be calculated nv <= N
rng(seed + 2, "twister");
duhat = zeros(Nx,nv);
dH = 0.1*rand(Nx,nv); % columns of dH are perturbation vectors
% dH = 0.001*orth(dH); % make them orthonormal and small
du = dH;
for iter = 1:nv
    duhat(:,iter) = fft(du(:,iter));
end
dH = du; % initializing perturbation vectors

laminst = zeros(nv,1); % lambda(k) instantaneous
dHmag = zeros(nv,nmax_ts); % magnitude of the perturbation vectors
dH1 = zeros(Nx,nmax_ts/nnorm);

tic
for n = 1:nmax_ts
    % dynamics----------:
    % -------predictor, corrector etd for u-------:
    uguess = utshat;
    for iters = 2
        upred = fft(real(ifft(utshat))) .* expalphadtts ...
                    + 0.5 * dt_tr * N1hat(uguess,kx) ...
                    + 0.5 * dt_tr * expalphadtts .* N1hat(utshat,kx);
        uguess = upred;
    end
    utshat = fft(real(ifft(utshat))) .* expalphadtts ...
                    + 0.5 * dt_tr * N1hat(upred,kx) ...
                    + 0.5 * dt_tr * expalphadtts .* N1hat(utshat,kx);
    
    u(:,n+1) = real(ifft(utshat));
    
    % tangent space dynamics-----------:
    % -------predictor, corrector etd for delta_u-------:
    % alphats = alpha - 1i*uts.*kx - ifft(1i*kx.*utshat);
    % for k = 1:nv
    %     duhat(:,k) = fft(real(ifft(duhat(:,k)))) .* exp(alphats*dt_ts);
    %     dH(:,k) = real(ifft(real(duhat(:,k))));
    % end
    % -------FOE for delta_u-------:
    alphats = alpha - 1i*uts.*kx - ifft(1i*kx.*utshat);
    for k = 1:nv
        duhat(:,k) = (alphats*dt_ts+1).*fft(real(ifft(duhat(:,k))));
        dH(:,k) = real(ifft(real(duhat(:,k))));
    end

    % renormalization and LLE ----------:
    if rem(n,nnorm) == 0
        [Q,R]= qr(dH);
        dH = Q(:,1:nv);
        dH1(:,n/nnorm) = Q(:,1);
        for k = 1:nv
            laminst(k,n/nnorm) = (1/tN)*log(abs(R(k,k)));
        end
    end
    
    % extracting dpsih and domzh from dH and calculating ||dH(k)||
    for k = 1:nv
        duhat(:,k) = fft(real(dH(:,k)));
        dHmag(k,n) = norm(dH(:,k));
    end

    du1(:,n) = dH(:,1);

    nnts = [nnts n];
    if rem(n,round(nmax_ts/10)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax_ts);
        toc;
    end
    if isnan(abs(sum(u)))
        error('blow up occured in du ..');
    end

end
toc

%%
figure; plot(dHmag(1,:),'-o');
figure; loglog(dHmag(1,:), '-o');
%%
figure; hold on;
for k = 1:nv
    plot(cumsum(laminst(k,:))./(1:length(laminst(k,:))), ...
        '-o', 'DisplayName', join(['\lambda' num2str(k)]));
end
legend;
% yline(mean(laminst(1,:)));
yline(0,'--');
%%
figure; 
subplot(1,2,1);
imagesc(u(:,nnorm:nnorm:end-1)'); colormap jet; colorbar;
set(gca,"YDir","normal");
subplot(1,2,2);
imagesc(abs(du1(:,nnorm:nnorm:end-1))');
clim([0 0.1])
colormap jet; colorbar;
set(gca,"YDir","normal");

%%
figure;
nnts = 1:nmax_ts;
surf(nnts(nnorm:nnorm:end-1),x,du1(:,nnorm:nnorm:end-1)), shading interp, lighting phong, axis tight 
view([-90 90]), colormap(autumn); set(gca,'zlim',[-5 5]) 
light('color',[1 1 0],'position',[-1,2,2]) 
material([0.30 0.60 0.60 40.00 1.00]);

%%
function N1h = N1hat(uhat,kx)
% function for the nonlinear part 
    N1h = -fft(real(ifft(uhat)).*real(ifft(1i*kx.*uhat)));
    % N1h = - 0.5*1i*kx.* fft(real(ifft(uhat)).^2 );
end

function N1hts = N1hts(uhat,duhat,kx)
N1hts = -fft(real(ifft(uhat)).*real(ifft(1i*kx.*duhat))) ...
        -fft(real(ifft(duhat)).*real(ifft(1i*kx.*uhat)));
end