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

% assigning v for fft of u
v = fft(uvec);

%% WAVENUMBER vector
kx = (1/spatial_res)*(-Nx/2:Nx/2-1)';
ky = (1/spatial_res)*(-Ny/2:Ny/2-1)';
kx = fftshift(kx);
ky = fftshift(ky);

[Kx,Ky] = meshgrid(kx,ky);

%% Linear derivative operator in the Fourier domain
L = -Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 + r - 1;


%% ETD explicit
tic
u = etd1exp_sh2D(u,v,dt,nmax,L,ToIntervals,Nx,Ny);

%%
figure; imagesc(u(:,:,end)); colormap jet; colorbar;
set(gca,'YDir','normal'); axis square;


%% wavenumber figure;
figure;
npad = 1000; % # of trailing zeros u is to be padded with
             % for fft calculation
f = spatial_res * ((-npad/2):(npad/2-1))/npad;
f = fftshift(f);
[X,Y] = meshgrid(f,f);
runsumYset = zeros(size(fft2(u(:,:,end),npad,npad)));
starttimestep = totimeu - round(3*totimeu/4);
endtimestep = totimeu;
for i = starttimestep:endtimestep
    Yset = fft2(u(:,:,i),npad,npad);
    runsumYset = runsumYset + Yset;
end
Yse = runsumYset./(endtimestep-starttimestep);
Yse = abs(Yse/sqrt(npad*npad)).^2;
contourf(X,Y,Yse,'LevelStep',1,'EdgeColor','none');
colorbar; colormap jet;
axis square;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
box("on");
xlabel('$k_x$','Interpreter','latex','FontSize',30);
ylabel('$k_y$','Interpreter','latex','FontSize',30);

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
