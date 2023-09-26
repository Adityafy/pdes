% 2D Swift Hohenberg
clear all;
close all;

clear all;
addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

% control parameter
eps = 0.7;
sig = 2;
c = sqrt(0.1);
gm = 50;
lam_0 = 2*pi; % critical wavelength

dt = 0.2; % time step size
dx = lam_0/8; % node spacing
dy = lam_0/8;

tmax = 1000/dt; % total time steps

rolls = 10;
Nx = round(rolls*lam_0/dx); % spatial nodes in x
Ny = round(rolls*lam_0/dy); % spatial nodes in y
N = Nx*Ny;
% N = 30;

seed = 2;
rng(seed,"twister");

para = [eps sig c gm dt dx dy tmax Nx Ny];

psilat = zeros(Nx,Ny,1); % scalar field
psivec = zeros(N,1);
psilat(:,:,1) = rand(Nx,Ny);
psivec(:,1) = latToVec(psilat(:,:,1));

zetalat = zeros(Nx,Ny,tmax); % scalar field
zetavec = zeros(N,tmax);
zetalat(:,:,1) = rand(Nx,Ny);
% zetalat(:,:,1) = zeros(Nx,Ny);
zetavec(:,1) = latToVec(zetalat(:,:,1));

omzlat = zeros(Nx,Ny,1); % scalar field
omzvec = zeros(N,1);
omzlat(:,:,1) = rand(Nx,Ny);
omzvec(:,1) = latToVec(zetalat(:,:,1));

tic

%% coefficients
z1 = 1/dx^2;
z2 = -2 * (1/dx^2 + 1/dy^2);
z3 = 1/dy^2;

a1 = dt*(eps-1)/2;
a2 = dt/(dx^4);
a3 = dt/(dx^2);
a4 = dt/(dx^2*dy^2);
a5 = dt/(dy^4);
a6 = dt/(dy^2);

a7 = a1 - 3*a2 - 4*a4 - 3*a5 + 2*a3 + 2*a6;

b1 = 1 - a7;
b2 = -2*a2 + a3 - 2*a4;
b3 = -2*a5 + a6 - 2*a4;
b4 = a4;
b5 = a2/2;
b6 = a5/2;

c1 = 1 + a7;
c2 = -b2;
c3 = -b3;
c4 = -b4;
c5 = -b5;
c6 = -b6;

m0 = sig*dt/(2*dx^2) + sig*dt/(2*dy^2) + c^2;
m1 = 1 + m0;
m2 = - sig*dt/(2*dx^2);
m3 = - sig*dt/(2*dy^2);
m4 = 1 - m0;
m5 = - m2;
m6 = - m3;

%% coefficient matrix construction
I = 1:Nx;
J = 1:Ny;
Ip1 = circshift(I,-1);
Im1 = circshift(I,1);
Ip2 = circshift(I,-2);
Im2 = circshift(I,2);
Jp1 = circshift(J,-1);
Jm1 = circshift(J,1);
Jp2 = circshift(J,-2);
Jm2 = circshift(J,2);

z = zeros(Nx,Ny);
Z = zeros(N,N);

b = zeros(Nx,Ny);
d = zeros(Nx,Ny);
B = zeros(N,N);
D = zeros(N,N);

ml = zeros(Nx,Ny);
mr = zeros(Nx,Ny);
Ml = zeros(N,N);
Mr = zeros(N,N);

n = 1;
for i = 1:length(I)
    for j = 1:length(J)
        z(I(i),J(j)) = z2;
        z(I(i),Jp1(j)) = z3;
        z(I(i),Jm1(j)) = z3;
        z(Ip1(i),J(j)) = z1;
        z(Im1(i),J(j)) = z1;
        Z(n,:) = latToVec(z)';

        b(I(i),J(j)) = b1;
        b(I(i),Jp1(j)) = b3;
        b(I(i),Jp2(j)) = b6;
        b(I(i),Jm1(j)) = b3;
        b(I(i),Jm2(j)) = b6;
        b(Ip1(i),J(j)) = b2;
        b(Ip1(i),Jp1(j)) = b4;
        b(Ip1(i),Jm1(j)) = b4;
        b(Ip2(i),J(j)) = b5;
        b(Im2(i),J(j)) = b5;
        b(Im1(i),J(j)) = b2;
        b(Im1(i),Jp1(j)) = b4;
        b(Im1(i),Jm1(j)) = b4;
        B(n,:) = latToVec(b)';
        d = -1*b;
        d(I(i),J(j)) = 1+a7;
        D(n,:) = latToVec(d)';
        
        ml(I(i),J(j)) = m1;
        ml(Im1(i),J(j)) = m2;
        ml(Ip1(i),J(j)) = m2;
        ml(I(i),Jm1(j)) = m3;
        ml(I(i),Jp1(j)) = m3;
        Ml(n,:) = latToVec(ml)';
        mr(I(i),J(j)) = 1 - m0;
        mr = - ml;
        Mr(n,:) = latToVec(mr)';

        n = n+1;
        z = zeros(Nx,Ny);
        b = zeros(Nx,Ny);
        d = zeros(Nx,Ny);
        ml = zeros(Nx,Ny);
        mr = zeros(Nx,Ny);
    end
end

Z = sparse(Z);

Lz = tril(Z);
Uz = triu(Z,1);
Lzinv = inv(Lz);

matdivpsi = B\D;
matdivomz = Ml\Mr;

%% Time integration
for t = 1:tmax
    zetavec(:,t) = Lzinv * (omzvec(:,t) - Uz * zetavec(:,t));
    zetalat(:,:,t) = vecToLat(zetavec(:,t),Nx,Ny);
    psitilde = rk2eq1(psilat(:,:,t),zetalat(:,:,t),para,@nlpgsh1); % explicit nonlinear
    psitilde = latToVec(psitilde);
    omztilde = rk2eq2(psilat(:,:,t),para,@nlpgsh2);
    omztilde = latToVec(omztilde);
    psivec(:,t+1) = matdivpsi*psitilde; % implicit CN linear
    omzvec(:,t+1) = matdivomz*omztilde;
    psilat(:,:,t+1) = vecToLat(psivec(:,t+1),Nx,Ny);
    omzlat(:,:,t+1) = vecToLat(omzvec(:,t+1),Nx,Ny);
    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
end
toc

%% Figures and videos

dynVideoFilename = 'gshDynVideo2';
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
%lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
figure();
for t = 1:50:tmax
    % imagesc(psilat(:,:,t));
    contourf(psilat(:,:,t), 'Linecolor', 'none');
    colorbar;
    colormap jet
    clim([-0.8 0.8]);
    %set(gca,'YDir','normal');
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
end
close(lat_dyn_video);


figure;
imagesc(psilat(:,:,tmax)');
set(gca,'YDir','normal');
colorbar;
colormap jet;

% figure;
% imagesc(ulat(:,:,6)');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;           

figure;
contourf(psilat(:,:,tmax)', 'Linecolor', 'none');
set(gca,'YDir','normal');
colorbar;
colormap jet;

figure;
imagesc(psivec');
set(gca,'YDir','normal');
colorbar; 
colormap jet;


%% Functions
function np1 = nlpgsh1(psi,zeta,para)
    dx = para(3);
    dy = para(4);
    Nx = para(9);
    Ny = para(10);
    %psi = vecToLat(psivec);
    %zeta = vecToLat(zetavec);
    I = 1:Nx;
    J = 1:Ny;
    Ip1 = circshift(I,-1);
    Im1 = circshift(I,1);
    Jp1 = circshift(J,-1);
    Jm1 = circshift(J,1);
    c = 1/(4*dx*dy);
    np1 = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            np1(i,j) = - 3*psi(I(i),J(j))^3 ...
                - c * (zeta(I(i),Jp1(j)) - zeta(I(i),Jm1(j))) ...
                    * (psi(Ip1(i),J(j)) - psi(Im1(i),J(j))) ...
                + c * (zeta(Ip1(i),J(j)) - zeta(Im1(i),J(j))) ...
                    * (psi(I(i),Jp1(j)) - psi(I(i),Jm1(j)));
        end
    end
end

function np2 = nlpgsh2(psi,para)
    gm = para(4);
    dx = para(6);
    dy = para(7);
    Nx = para(9);
    Ny = para(10);
    %psi = vecToLat(psivec);
    % zeta = vecToLat(zetanvec);
    I = 1:Nx;
    J = 1:Ny;
    Ip1 = circshift(I,-1);
    Im1 = circshift(I,1);
    Ip2 = circshift(I,-2);
    Im2 = circshift(I,2);
    Jp1 = circshift(J,-1);
    Jm1 = circshift(J,1);
    Jp2 = circshift(J,-2);
    Jm2 = circshift(J,2);
    
    dpsidx = @(i,j) (psi(Ip1(i),J(j)) - psi(Im1(i),J(j)))/(2*dx);
    %dzetadx = @(i,j) (zeta(Ip1(i),J(j)) - zeta(Im1(i),J(j)))/(2*dx);
    dpsidy = @(i,j) (psi(Ip1(i),J(j)) - psi(Im1(i),J(j)))/(2*dy);
    %dzetady = @(i,j) (zeta(Ip1(i),J(j)) - zeta(Im1(i),J(j)))/(2*dy);
    
    d3psidx3 = @(i,j) (1/2*dx^3) * (-psi(Im2(i),J(j)) + 2 * psi(Im1(i),J(j)) ...
            - 2 * psi(Ip1(i),J(j)) + psi(Ip2(i),J(j)));
    d3psidy3 = @(i,j) (1/2*dy^3) * (-psi(I(i),Jm2(j)) + 2 * psi(I(i),Jm1(j)) ...
            - 2 * psi(I(i),Jp1(j)) + psi(I(i),Jp2(j)));
    d3psidxdy2 = @(i,j) (1/dx*dy^2) * (psi(Ip1(i),Jm1(j)) - psi(Im1(i),Jm1(j)) ...
            - 2 * psi(Ip1(i),J(j)) + 2*psi(Im1(i),J(j)) ...
            + psi(Ip1(i),Jp1(j)) - psi(Im1(i),Jp1(j)));
    d3psidydx2 = @(i,j) (1/dy*dx^2) * (psi(Im1(i),Jp1(j)) - psi(Im1(i),Jm1(j)) ...
            - 2 * psi(Ip1(i),J(j)) + 2*psi(Im1(i),J(j)) ...
            + psi(Ip1(i),Jp1(j)) - psi(Im1(i),Jp1(j)));
    np2 = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            np2(i,j) = gm * ( dpsidy(i,j) * (d3psidx3(i,j) + d3psidxdy2(i,j)) ...
                    - dpsidx(i,j) * ( d3psidydx2(i,j) + d3psidy3(i,j) ) );
        end
    end
end

function psi = rk2eq1(psi,zeta,para,dynfunc)
    dt = para(5);
    k1 = dynfunc(psi,zeta,para);
    %u1 = u+k1*dt;
    k2 = dynfunc(psi+k1*dt,zeta,para);
    psi = psi + dt*((k1+k2)/2);
end

function psi = rk2eq2(psi,para,dynfunc)
    dt = para(5);
    k1 = dynfunc(psi,para);
    %u1 = u+k1*dt;
    k2 = dynfunc(psi+k1*dt,para);
    psi = psi + dt*((k1+k2)/2);
end

function u = rk4(u,dynFunc,h)
    % k values
    k1 = dynFunc(u);
    k2 = dynFunc(u + (0.5*h)*k1);
    k3 = dynFunc(u + (0.5*h)*k2);
    k4 = dynFunc(u + h*k3);
    % dynamics
    u = u + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
end

function xvec = latToVec(x)
    [Nx,Ny] = size(x);
    xvec = zeros(Nx*Ny,1);
    for i = 1:Nx
        for j = 1:Ny
            k = (i-1)*Ny+j;
            xvec(k) = x(i,j);
        end
    end
    % xvec = xvec';
    % xdiag = diag(xvec);
end

function xlat = vecToLat(x,Nx,Ny)
    if length(x) ~= Nx*Ny
        msg = "Output lattice dimensions do not comply " + ...
            "with the length of input vector";
        error(msg);
    end
    xlat = zeros(Nx,Ny);
    for i = 1:Nx*Ny
        q = fix(i/Ny);
        r = rem(i,Ny);
        k1 = q+1;
        k2 = r;
        if r == 0
            k1 = q;
            k2 = Ny;
        end
        % xlat(k1,k2) = x(i,i);
        xlat(k1,k2) = x(i);
    end
end
