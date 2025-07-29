%% Testing GMRES and other iterative methods
% First we are going to test laplacian operator itself.
% -nabla^2 zeta = Omega_z
% We initiate zeta first and calculate the Omega_z using finite difference
% and spectral methods.
% After that we calculate zeta from various iterative methods using
% omz and compares it with the original initiated value.
clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% spatial parameters
dsm = 2; % domain size multiple
resoln = 8; % resolution should be >1. 8 and more is very good.
Lx = 2*pi*dsm; % Integral multiple of 2*pi for functions with period 2*pi.
Ly = Lx; % Better keep it this way
L = Lx; % domain length
Nx = ceil(resoln*Lx); % this way Nx is always more than L
if rem(Nx,2)~=0
    Nx = Nx+1; % this way Nx is always an even #
end
Ny = Nx;
N = Nx*Ny; % system size

% ---------  spatial grid in real domain -----------
x2 = linspace(-Lx/2,Lx/2,Nx+1)';
x = x2(1:Nx);
y2 = linspace(-Ly/2,Ly/2,Ny+1)';
y = y2(1:Ny);
dx = x(2) - x(1);
dy = y(2) - y(1);

[X,Y] = meshgrid(x,y);

%% Wavenumbers

kx = (2*pi/Lx)*(-Nx/2:Nx/2-1)';
ky = (2*pi/Ly)*(-Ny/2:Ny/2-1)';
kx = fftshift(kx);
ky = fftshift(ky);
% kx(1) = 10^-2;
% ky(1) = 10^-2;
% eps = 0;
[Kx,Ky] = meshgrid(kx,ky);
% Kx = Kx';
% Ky = Ky';

Kdiff = Kx.^2 + Ky.^2;
firstrowKdiff = latToVec(Kdiff)';
KdiffMat = toeplitz(firstrowKdiff);
qx = latToVec(Kx);
qy = latToVec(Ky);

%% zeta and omz

zetamat = (cos(2*X)-sin(2*Y));
delzeta = 2*(-sin(2*X)-cos(2*Y));
omzmattheory = -4*(-cos(2*X)+sin(2*Y));

% zetamat = (cos(X)+sin(Y));
% delzeta = (-sin(X)+cos(Y));
% omzmattheory = -(-cos(X)-sin(Y));

% zetamat = cos(X/2);
% delzeta = -0.5*sin(X/2);
% omzmattheory = -(-0.25*cos(X/2));

%%
figure;
subplot(3,1,1); plot(x,zetamat(1,:),'-o');
subplot(3,1,2); plot(x,delzeta(1,:),'-o');
subplot(3,1,3); plot(x,omzmattheory(1,:),'-o');

%%
figure; imagesc(x,y,zetamat); title('\zeta initial'); colorbar;
set(gca,'YDir','normal');
%%
figure; imagesc(x,y,omzmattheory); title('\Omega theory'); colorbar;
set(gca,'YDir','normal');
%% finite difference omz
omzmatfd = -fdlaplacian(zetamat,Nx,Ny,dx);
%%
figure; imagesc(x,y,omzmatfd); title('\Omega finite-difference'); 
colorbar; set(gca,'YDir','normal');

%%
zeta = latToVec(zetamat);
% umat = zeros(length(X), length(Y),1);
% vmat = zeros(length(X), length(Y),1);
% utr = zeros(length(X), length(Y),1);
% vtr = zeros(length(X), length(Y),1);
zetavec = fftshift(zeta);
zetah = fft(zetavec);
zeta2h = fft2(zetamat);

%%
% omztr2h = -Kdiff.*fft2(zetatrmat);
% omztrmat = real(ifft2(omztr2h));
omz2h = Kdiff.*(zeta2h);
omzmatspec = real(ifft2(omz2h));
omzh = fft(latToVec(omzmatspec));
omzspec = real(ifft(omzh));

%%
% figure; imagesc(kx,ky,real(omz2h)); colorbar;
%%
figure; imagesc(x,y,omzmatspec); title('\Omega spectral'); colorbar;
set(gca,'YDir','normal');

%%
npad = 1000; % # of trailing zeros psi is to be padded with
             % for fft calculation
f = (2*pi/Lx)*((-npad/2):(npad/2))/(npad+1);
f = f(1:end-1);
f = fftshift(f);
[Xpsi,Ypsi] = meshgrid(f,f);
Yset = fft2(zetamat,npad,npad);
Yset = abs(Yset/sqrt(npad*npad)).^2;

figure;
% contourf(Xpsi,Ypsi,Yset,'LevelStep',0.1,'EdgeColor','none');
contourf(Xpsi,Ypsi,Yset);
colorbar; colormap jet;
axis square;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
box("on");
xlabel('$k_x$','Interpreter','latex','FontSize',30);
ylabel('$k_y$','Interpreter','latex','FontSize',30);

%% zeta fd
zeta_fd = myIterZetaOmz(omzmattheory,zeros(size(zetamat)),Nx,Ny,dx);

figure; imagesc(zeta_fd); colorbar; axis square; 
title('\zeta finite difference');

%% zeta gmres
% Ax = b => fdA * zetavec_gmres = omzvecfd

%% In built solvers

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

%%% converting omz calculated by fd to vector form
omzvecfd = latToVec(omzmatfd);

%%% calling gmres
restart = 100;
tol = 1e-12; % error tolerance
maxit = 500; % max # of iterations
zetavec_gmres = gmres(fdA,omzvecfd,restart,tol,maxit);
% zetavec_gmres = qmr(fdA',omzvecfd,tol,maxit);

zeta_gmres = vecToLat(zetavec_gmres,Nx,Ny);

%%
figure;
subplot(1,3,1); imagesc(x,y,zetamat); title('\zeta initial');
colorbar; axis square;
subplot(1,3,3); imagesc(y,x,zeta_gmres); title('\zeta gmres');
colorbar; axis square;

%%
figure; imagesc(zetamat - zeta_gmres); colorbar;
%%
function g = fdlaplacian(f, Nx, Ny, dx)
% g = fdlaplacian(f, Nx, Ny, dx)
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

%%
function zeta = myIterZetaOmz(omzmat, zetamat, Nx,Ny,dx)
% zeta = myIterZetaOmz(omzmat, zetamat, Nx,Ny,dx)
% An iterative poisson solver for
% Omega_z = - Laplacian(Zeta)
% with periodic boundary conditions.
% mind the negative sign if used for other kinds of Poisson solvers
% Number of iterations is fixed (change it in function)
    delta = dx;

    I = 1:Nx;
    J = 1:Ny;
    Ip1 = circshift(I,-1);
    Im1 = circshift(I,1);
    Jp1 = circshift(J,-1);
    Jm1 = circshift(J,1);

    zeta = zetamat;
    iterations = 800;
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
