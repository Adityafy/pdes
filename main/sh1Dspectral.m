%% Swift-Hohenberg 1D with Pseudospectral and ETD
% Now works with different L (domain) and N (grid points)
% using various time-integration solvers:
% First order Euler. ETD1 explicit (predictor corrector).
% SHE in 1D is: u_t = (r-1)u - 2u_xx - u_xxxx - u^3.

clear all;
close all;

addpath('../src/');

mfa = '../../pdesDataDump/media/'; % media folder address
dfa = '../../pdesDataDump/data/'; % saved data folder address

%% Parameters
r = 0.05; % control parameter
lam_0 = 2*pi; % critical wavelength

%% SPACE
Lx = 12; % length of the square domain
Nx = 96; % number of grid points

% ---------  spatial grid in real domain -----------
x2 = linspace(-Lx/2,Lx/2,Nx+1)';
x = x2(1:Nx);
dx = x(2) - x(1);

%% TIME
dt = 0.2; % time step size
totimeu = 500;
nmax = totimeu/dt; % total time steps

%% INITIAL CONDITIONS

% -----periodic ICs-----
% u = cos(x/16).*(1+sin(x/16)); 

% -----random ICs-----
seed = 2;
rng(seed,"twister");
u = zeros(size(x));
u(:,1) = 0.01*rand(size(u));

%% Analytical solution
A1 = sqrt(4*r/3);
% A3 = -A1/4 + sqrt(A1^2/16 + r*A1/3 + -A1^3/4);
A3 = r^(3/2);
A5 = r^(5/2);
u_solution = A1*cos(x) + A3*cos(3*x);


%% WAVENUMBER vector
kx = (2*pi/Lx)*(-Nx/2:Nx/2-1)'; % wavenumber grid
kx = fftshift(kx);
L = r - 1 + 2*kx.^2 - kx.^4; % Fourier multipliers

%% First order Euler
% tic
% for n = 1:nmax
%     k1 = L.*fft(real(ifft(v)))+fft(real(-(ifft(v)).^3));
%     v = v + dt * k1;
%     u(:,n+1) = real(ifft(v));
%     if rem(n,round(nmax/10)) == 0
%         % clc;
%         fprintf('This is time step: %g / %g', n, nmax);
%         toc;
%     end
% end
% toc


% % Space time dynamics figure
% figure;
% % contourf(u','LevelStep',0.1,'EdgeColor','none'); 
% imagesc(u');
% colormap jet; colorbar;
% set(gca,'YDir','normal');


%% ETD explicit
v = fft(u(:,1));
tic
uetd1exp = etd1exp_sh1D(v,dt,nmax,L);

figure;
% contourf(u','LevelStep',0.1,'EdgeColor','none'); 
imagesc(uetd1exp');
colormap jet; colorbar;
set(gca,'YDir','normal');

%% 
figure; plot(uetd1exp(1,:),'-o','LineWidth',1);
xlabel('$n$','FontSize',30,'Interpreter','latex');
ylabel('$u^{(1)}$','FontSize',30,'Interpreter','latex','Rotation',0);

%% wavenumber
figure;
ufftreal = real(fft(uetd1exp(:,end)));
plot(kx,ufftreal,'-o','LineWidth',1);
xlim([min(kx) max(kx)]);
set(gca,'TickLabelInterpreter','tex','FontSize',15);
axis square;
box("on");
xlabel('$k$','Interpreter','latex','FontSize',30);
ylabel('$\hat{u}_{ss}$','Interpreter','latex','FontSize',30,'Rotation',0);

%% functions

%% ETD1 function
function u = etd1exp_sh1D(v,dt,nmax,L)
% This is done by the predictor corrector method used in Cross et al.
% (1994). There is an extra loop in the time integration loop so that
% predictor guess is converged.
u = zeros(size(real(ifft(v))));
for n = 1:nmax
    % predictor1
    v_pred = v; % initial guess for predictor step
    % v_pred = fft(real(ifft(rand(size(u(:,1))))));
    for i = 1
        v_pred = fft(real(ifft(v))).*exp(L*dt) ...
                    + 0.5*dt*fft(real(-(ifft(v_pred)).^3)) ...
                       + 0.5*dt*exp(L*dt).*fft(real(-(ifft(v)).^3));
    end
    % corrector (final step for u_{n+1})
    v = fft(real(ifft(v))).*exp(L*dt) + 0.5*dt*fft(real(-(ifft(v_pred)).^3)) + ...
                        0.5*dt*exp(L*dt).*fft(real(-(ifft(v)).^3));
    u(:,n+1) = real(ifft(v));
    if rem(n,round(nmax/10)) == 0
        
        fprintf('This is time step: %g / %g, ', n, nmax);
        toc;
    end
end

end

%% ETD RK4 function
function u = etdrk4(p,u,Lfunc,Nfunc)
% u = etdCoxMatthewsRK4(p,u)

h = p.dt;
m = p.m;
nmax = p.nmax;
t = p.time;

% precomputing some
L = Lfunc(p);
N = Nfunc;

E1 = expm(L*h);
E2 = expm(L*h/2);
% E1 = exp(L*h);
% E2 = exp(L*h/2);
% omega = inv(L) * (E2-eye(m));
omega = L \ (E2-eye(m));

eta = h^(-2) * L^(-3);
alpha = eta * (- 4 - L*h + E1 * (4 - 3*L*h + (L*h)^2));
beta = eta * (2 + L*h + E1 * (-2+L*h));
gamma = eta * (-4 - 3*L*h - (L*h)^2 + E1 * (4 - L*h));

% time integration
for n = 1:nmax
    an = E2 * u(:,n) + omega * N(u(:,n),t(n));
    bn = E2 * u(:,n) + omega * N(an,t(n)+h/2);
    cn = E2 * an + omega * ( 2 * N(bn,t(n)+h/2) - N(u(:,n),t(n)) );

    u(:,n+1) =  E1 * u(:,n) + alpha * N(u(:,n), t(n)) ...
        + beta * 2 * ( N(an, t(n)+h/2) + N(bn, t(n)+h/2) ) ...
        + gamma * N(cn, t(n)+h) ;

end

end