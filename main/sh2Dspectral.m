%% Script for Swift-Hohenberg 2D with Spectral
% using various time-integration solvers. Particularly ETD.
% First order Euler works.
% ETD1 explicit (predictor corrector) also works.
% SHE in 2D is: 
%   u_t = (r-1)u - u_xxxx - 2u_xxyy - u_yyyy - 2u_xx - 2u_yy - u^3.

clear all;
close all;

addpath('../src/');
mfa = '../../pdesDataDump/media/'; % media folder address
dfa = '../../pdesDataDump/data/'; % saved data folder address

%% Parameters
r = 0.7; % control parameter
lam_0 = 2*pi; % critical wavelength

%% SPACE
spatial_res = 16;

dx = lam_0/spatial_res; % node spacing

Nx = spatial_res^2;
Ny = Nx;

x = spatial_res*2*pi*(1:Nx)'/Nx;
y = spatial_res*2*pi*(1:Ny)'/Ny;

% making the grid
[X,Y] = meshgrid(x,y);

%% TIME
dt = 0.2; % time step size
totimeu = 50;
nmax = totimeu/dt; % total time steps
ToIntervals = totimeu;

%% INITIAL CONDITIONS

% -----periodic ICs-----
% u = cos(x/spatial_res).*(1+sin(x/spatial_res)); 

% -----random ICs-----
seed = 2;
rng(seed,"twister");

u(:,:,1) = 0.1*rand(length(X), length(Y));

% changing to vector
uvec = latToVec(u);

% assigning v for fft of u
v = fft(uvec);

%% WAVENUMBER vector
kx = [0:Nx/2-1 0 -Nx/2+1:-1]'/spatial_res; % wave numbers in x
ky = [0:Ny/2-1 0 -Ny/2+1:-1]'/spatial_res; % wave numbers in y
% kx = [0:(Nx/2-1) (-Nx/2):-1]'/spatial_res; % wave numbers in x
% kx(1) = 10^(-6);
% ky = [0:(Ny/2-1) (-Ny/2):-1]'/spatial_res; % wave numbers in y
% ky(1) = 10^(-6);
[Kx,Ky] = meshgrid(kx,ky);

% the linear operator matrix after the fourier transform
L = -Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 + r - 1;


%% ETD explicit
% v = fft(u(:,1));
tic
u = etd1exp_sh2D(u,v,dt,nmax,L,ToIntervals,Nx,Ny);

%%
figure; imagesc(u(:,:,end)); colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;

%% wavenumber
% figure; hold on;
% plot(latToVec(Kx),real(fft(latToVec(u(:,:,end)))),'-o','LineWidth',1);
% % xlim([min(k_vec) max(k_vec)]);
% set(gca,'TickLabelInterpreter','tex','FontSize',15);
% axis square;
% box("on");
% xlabel('$k$','Interpreter','latex','FontSize',30);
% ylabel('$\hat{u}_{ss}$','Interpreter','latex','FontSize',30);

%% functions

%% ETD1 function
function u = etd1exp_sh2D(u,v,dt,nmax,L,ToIntervals,Nx,Ny)
% This is done by the predictor corrector method used in Cross et al.
% (1994). There is an extra loop in the time integration loop so that
% predictor guess is converged.

expLhvec = exp(latToVec(L*dt));
for n = 1:nmax
    % -------------predictor-------------
    v_pred = v; % initial guess for predictor step
    % v_pred = fft(real(ifft(rand(size(u(:,1))))));
    for i = 1
        v_pred = fft(real(ifft(v))).*expLhvec ...
                    + 0.5*dt*fft(real(-(ifft(v_pred)).^3)) ...
                       + 0.5*dt*expLhvec.*fft(real(-(ifft(v)).^3));
    end
    % -------------corrector (final step for u_{n+1})-------------
    v = fft(real(ifft(v))).*expLhvec + 0.5*dt*fft(real(-(ifft(v_pred)).^3)) + ...
                        0.5*dt*expLhvec.*fft(real(-(ifft(v)).^3));

    if rem(n,nmax/ToIntervals) == 0
        %------------- saving dynamics at intervals ---------------
        umat = vecToLat(real(ifft(v)),Nx,Ny);
        u(:,:,round(n*ToIntervals/nmax)) = umat;
    end
    if rem(n,round(nmax/10)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax);
        toc;
    end
end

end
