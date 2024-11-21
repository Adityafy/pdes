% 1D Swift Hohenberg

close all;

clear all;
addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

r = 0.1; % control parameter

lam_0 = 2*pi; % critical wavelength

dt = 0.01; % time step size
% dx = 0.99;
spatial_res = 16;
dx = lam_0/spatial_res; % node spacing


totimeu = 2000;
nmax = totimeu/dt; % total time steps

% rolls = 32;
% N = round(rolls*lam_0/dx); % spatial nodes

N = 256;
k_vec = [0:N/2-1 0 -N/2+1:-1]'/spatial_res; % wave numbers 

%%
seed = 2;
rng(seed,"twister");

para = [r dt dx nmax N];

u = zeros(N,1); % scalar field
u(:,1) = 0.01*rand(N,1);
% u(:,1) = zeros(N,1);
% u(:,1) = ones(N,1);

tic

%% coefficient matrix construction for the linear part L[u]
c1 = dt/dx^2;
c2 = dt/dx^4;

ids = 1; % implicit division steps

% LHS
a2 = c1 - 2*c2;
a4 = a2;
a3 = ids - dt*(1/2)*(r-1) - 2*c1 + 3*c2;
a1 = c2/2;
a5 = a1;
fl = c2;

cfmat_lhs = zeros(N,N);
cvec_lhs = zeros(1,N);
cvec_lhs(1,1:3) = [a3 a4 a5];
cvec_lhs(1,end-1:end) = [a1 a2];
cfmat_lhs(1,:) = cvec_lhs;

% RHS

b3 = ids + dt*(1/2)*(r-1) + 2*c1 - 3*c2;
b2 = - a2;
b4 = b2;
b1 = - a1;
b5 = b1;
br = - c2;

cfmat_rhs = zeros(N,N);
cvec_rhs = zeros(1,N);
cvec_rhs(1,1:3) = [b3 b4 b5];
cvec_rhs(1,end-1:end) = [b1 b2];
cfmat_rhs(1,:) = cvec_rhs;

for i = 2:N
    cfmat_rhs(i,:) = circshift(cvec_rhs,i-1);
    cfmat_lhs(i,:) = circshift(cvec_lhs,i-1);
end

matdiv = inv(cfmat_lhs)*cfmat_rhs;

%% Time integration
for npad = 1:nmax
    ustar = rk2(u(:,npad),@nonlinPart,dt); % explicit nonlinear
    %ustar2 = matdiv*ustar;
    u(:,npad+1) = matdiv*ustar; % implicit cn linear
    if rem(npad,nmax/10) == 0
        fprintf('This is time step %i, ',npad);
        toc
    end
end
toc

%%
figure;
imagesc(u');
set(gca,'YDir','normal');
colorbar;
colormap jet;


%%
% x = linspace(0,dx*(N),N);
% u_solution = sqrt(4*r/3)*cos(x) + (((sqrt(4/3))^3 * cos(3*x))./(4*(r-64)));
x = 0:dx:dx*N;

A1 = sqrt(4*r/3);
% A3 = -A1/4 + sqrt(A1^2/16 + r*A1/3 + -A1^3/4);
A3 = 0.14*A1*r;
A5 = 0.05*A3*r;
u_solution = A1*cos(x) + A3*cos(3*x) + A5*cos(5*x);
% N = 30;
figure;
plot(u(:,end),'-o', 'DisplayName','finite difference solution');
hold on;
plot(u_solution,'DisplayName','analytical'); legend;

%%
% figure;
% contourf(u');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;
% 
% figure;
% imagesc(u(:,1:50)');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;
% 
% figure;
% imagesc(u(:,end-500:end)');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;

%% Is it steady state?
figure;
plot(u(25,:),'-o','DisplayName','u_1(t)');
xlabel('t');
ylabel('u_1');
subtitle('Is the solution steady state?');
xlim([1 length(u(1,:))]);

%% wavenumber
figure; hold on;

npad = 10000; % # of trailing zeros u is to be padded with 
             % for fft calculation
Y = fft(u(:,end),npad);
P = abs(Y/sqrt(npad)).^2;
f = spatial_res * (0:(npad/2))/npad;
plot(f,P(1:npad/2+1),'-o');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
axis square;
box("on");
xlabel('$k$','Interpreter','latex','FontSize',30);
ylabel('$|\hat{u}|$','Interpreter','latex','FontSize',30);

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

