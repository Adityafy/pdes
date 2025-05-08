%% Test Script for Swift-Hohenberg 2D with Pseudospectral
% Now works with different L (domain) and N (grid points)
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

Lx = 96;
Ly = Lx;
Nx = 256;
Ny = Nx;

% ---------  spatial grid in real domain -----------
x2 = linspace(-Lx/2,Lx/2,Nx+1)';
x = x2(1:Nx);
y2 = linspace(-Ly/2,Ly/2,Ny+1)';
y = y2(1:Ny);
dx = x(2) - x(1);
dy = y(2) - y(1);

[X,Y] = meshgrid(x,y);

%% TIME
dt = 0.1; % time step size
totimeu = 100;
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
umat = u;

% assigning v for fft of u
% v = fft(uvec);
vmat = fft2(umat);

%% wavenumber vector

kx = (2*pi/Lx)*(-Nx/2:Nx/2-1)';
ky = (2*pi/Ly)*(-Ny/2:Ny/2-1)';
kx = fftshift(kx);
ky = fftshift(ky);
[Kx,Ky] = meshgrid(kx,ky);

% the linear operator matrix after the fourier transform
L = -Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 + r - 1;


%% ETD explicit
% v = fft(u(:,1));
tic
% u = etd1exp_sh2D(u,v,dt,nmax,L,ToIntervals,Nx,Ny);
expLhmat = vecToLat(exp(latToVec(L*dt)),Nx,Ny);
for n = 1:nmax
    % -------------predictor-------------
    v_pred = vmat; % initial guess for predictor step
    % v_pred = fft(real(ifft(rand(size(u(:,1))))));
    
    for i = 1
        v_pred = fft2(real(ifft2(vmat))).*expLhmat ...
                    + 0.5*dt*fft2(real(-(ifft2(v_pred)).^3)) ...
                       + 0.5*dt*expLhmat.*fft2(real(-(ifft2(vmat)).^3));
    end
    % -------------corrector (final step for u_{n+1})-------------
    vmat = fft2(real(ifft2(vmat))).*expLhmat ...
            + 0.5*dt*fft2(real(-(ifft2(v_pred)).^3)) ...
                      + 0.5*dt*expLhmat.*fft2(real(-(ifft2(vmat)).^3));

    if rem(n,nmax/ToIntervals) == 0
        %------------- saving dynamics at intervals ---------------
        umat = real(ifft2(vmat));
        u(:,:,round(n*ToIntervals/nmax)) = umat;
    end
    if rem(n,round(nmax/10)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax);
        toc;
    end
end

%%
figure; imagesc(u(:,:,end)); colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;

%% wavenumber
figure; hold on;
plot(latToVec(Kx),real(fft(latToVec(u(:,:,end)))),'-o','LineWidth',1);
% xlim([min(k_vec) max(k_vec)]);
set(gca,'TickLabelInterpreter','tex','FontSize',15);
axis square;
box("on");
xlabel('$k$','Interpreter','latex','FontSize',30);
ylabel('$\hat{u}_{ss}$','Interpreter','latex','FontSize',30);

%% functions

%% ETD1 function
% function u = etd1exp_sh2D(u,v,dt,nmax,L,ToIntervals,Nx,Ny)
% % This is done by the predictor corrector method used in Cross et al.
% % (1994). There is an extra loop in the time integration loop so that
% % predictor guess is converged.
% 
% expLhvec = exp(latToVec(L*dt));
% for n = 1:nmax
%     % -------------predictor-------------
%     v_pred = v; % initial guess for predictor step
%     % v_pred = fft(real(ifft(rand(size(u(:,1))))));
%     for i = 1
%         v_pred = fft(real(ifft(v))).*expLhvec ...
%                     + 0.5*dt*fft(real(-(ifft(v_pred)).^3)) ...
%                        + 0.5*dt*expLhvec.*fft(real(-(ifft(v)).^3));
%     end
%     % -------------corrector (final step for u_{n+1})-------------
%     v = fft(real(ifft(v))).*expLhvec + 0.5*dt*fft(real(-(ifft(v_pred)).^3)) + ...
%                         0.5*dt*expLhvec.*fft(real(-(ifft(v)).^3));
% 
%     if rem(n,nmax/ToIntervals) == 0
%         %------------- saving dynamics at intervals ---------------
%         umat = vecToLat(real(ifft(v)),Nx,Ny);
%         u(:,:,round(n*ToIntervals/nmax)) = umat;
%     end
%     if rem(n,round(nmax/10)) == 0
%         fprintf('This is time step: %g / %g, ', n, nmax);
%         toc;
%     end
% end
% 
% end
