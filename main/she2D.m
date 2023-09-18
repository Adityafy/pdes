% 2D Swift Hohenberg

close all;

clear all;
addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

r = 0.3; % control parameter

lam_0 = 2*pi; % critical wavelength

dt = 0.2; % time step size
dx = lam_0/8; % node spacing
dy = lam_0/8;

tmax = 10^2/dt; % total time steps

rolls = 13;
Nx = round(rolls*lam_0/dx); % spatial nodes in x
Ny = round(rolls*lam_0/dy); % spatial nodes in y
N = Nx*Ny;
% N = 30;

seed = 2;
rng(seed,"twister");

para = [r dt dx dy tmax Nx Ny];

ulat = zeros(Nx,Ny,1); % scalar field
uvec = zeros(N,1);
ulat(:,:,1) = rand(Nx,Ny);
uvec(:,1) = latToVec(ulat(:,:,1));

tic

%% coefficients
a1 = dt*(r-1)/2;
a2 = dt/(dx^4);
a3 = dt/(dx^2);
a4 = dt/(dx^2*dy^2);
a5 = dt/(dy^4);
a6 = dt/(dy^2);

%a8 = -a1 - 3*a2 -4*a4 + 3*
% l= a1 - 3*a2 - 4*a4 - 3*a5 + 2*a3 + 2*a6;
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

b = zeros(Nx,Ny);
c = zeros(Nx,Ny);
B = zeros(N,N);
C = zeros(N,N);

n = 1;
for i = 1:length(I)
    for j = 1:length(J)
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
        c = -1*b;
        c(I(i),J(j)) = 1+a7;
        C(n,:) = latToVec(c)';
        n = n+1;
        b = zeros(Nx,Ny);
        c = zeros(Nx,Ny);
    end
end


% 
% b(1,1) = b1;
% b(1,2) = b3;
% b(1,3) = b6;
% b(1,Ny) = b2;
% b(1,Ny-1) = b6;
% b(2,1) = b3;
% b(2,2) = b4;
% b(2,end) = b4;
% b(3,1) = b5;
% b(Nx-1,1) = b5;
% b(Nx,1) = b2;
% b(Nx,2) = b4;
% b(Nx,Ny) = b4;
% 
% c = -1*b;
% c(1,1) = 1 + a7;
% 
% bvec = latToVec(b)';
% cvec = latToVec(c)';
% 
% 
% 
% for i = 1:N
%     B(i,:) = circshift(bvec,i-1);
%     C(i,:) = circshift(cvec,i-1);
% end

matdiv = B\C;

%% Time integration
for t = 1:tmax
    utilde = rk2(uvec(:,t),@nonlinPart,dt); % explicit nonlinear
    uvec(:,t+1) = matdiv*utilde; % implicit CN linear
    ulat(:,:,t+1) = vecToLat(uvec(:,t+1),Nx,Ny);
    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
end
toc



figure;
imagesc(ulat(:,:,tmax)');
set(gca,'YDir','normal');
colorbar;
colormap jet;

% figure;
% imagesc(ulat(:,:,6)');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;           

figure;
contourf(ulat(:,:,tmax)', 'Linecolor', 'none');
set(gca,'YDir','normal');
colorbar;
colormap jet;

figure;
imagesc(uvec');
set(gca,'YDir','normal');
colorbar; 
colormap jet;


%% Functions
function Nu = nonlinPart(u)
    Nu = - u.^3;
end

function u = rk2(u,dynfunc,dt)
    k1 = dynfunc(u);
    %u1 = u+k1*dt;
    k2 = dynfunc(u+k1*dt);
    u = u + dt*((k1+k2)/2);
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

