%% Test script for Swift-Hohenberg 1D with Spectral
% using various time-integration solvers. Particularly ETD.
% First order Euler works.
% ETD1 explicit (predictor corrector) also works.
% SHE in 1D is: u_t = (r-1)u - 2u_xx - u_xxxx - u^3.

clear all;
addpath('../src/');

mfa = '../../pdesDataDump/media/'; % media folder address
dfa = '../../pdesDataDump/data/'; % saved data folder address

%% Parameters
r = 0.1; % control parameter
lam_0 = 2*pi; % critical wavelength

%% SPACE
% I saw that spatial resolution of 16 works best in the follwing setting.
spatial_res = 16;

dx = lam_0/spatial_res; % node spacing
% num_rolls = 16;
% N = round(num_rolls*spatial_res);

N = spatial_res^2;
% x = 32*2*pi*(1:N)'/N;

x = spatial_res*2*pi*(1:N)'/N;
% x = spatial_res*2*pi*([-N/2:-1 0 1:N/2-1])'/N;

%% TIME
dt = 0.2; % time step size
totimeu = 2000;
nmax = totimeu/dt; % total time steps

%% INITIAL CONDITIONS

% -----periodic ICs-----
% u = cos(x/16).*(1+sin(x/16)); 

% -----random ICs-----
seed = 2;
rng(seed,"twister");
u = zeros(size(x));
u(:,1) = 0.01*rand(size(u));

% -----In Fourier space----
v = fft(u(:,1)); % initializing Fourier space variable


%% Analytical solution
x = linspace(0,dx*(N),N);
% u_solution = sqrt(4*r/3)*cos(x) + (((sqrt(4/3))^3 * cos(3*x))./(4*(r-64)));
A1 = sqrt(4*r/3);
% A3 = -A1/4 + sqrt(A1^2/16 + r*A1/3 + -A1^3/4);
A3 = r^(3/2);
A5 = r^(5/2);
u_solution = A1*cos(x) + A3*cos(3*x);
% N = 30;

%% WAVENUMBER vector
% k_vec = [0:N/2-1 0 -N/2+1:-1]'/32; % wave numbers 
% k_vec = [0:N/2-1 0 -N/2+1:-1]'/spatial_res; % wave numbers 
k_vec = (1/spatial_res)*(-N/2:N/2-1)';
k_vec = fftshift(k_vec);
% k_vec = [0:N/2-1 0 -N/2+1:-1]'; % wave numbers 
% k_vec = [-N/2:-1 0 1:N/2-1]'/spatial_res; % wave numbers
L = r - 1 + 2*k_vec.^2 - k_vec.^4; % Fourier multipliers

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


%% Space time dynamics figure
% figure;
% % contourf(u','LevelStep',0.1,'EdgeColor','none'); 
% imagesc(u');
% colormap jet; colorbar;
% set(gca,'YDir','normal');


%% ETD explicit
tic
uetd1exp = etd1exp_sh1D(v,dt,nmax,L);

%% figures
% figure;
% plot(u(:,end),'-','DisplayName','FOE');
% hold on;
% plot(uetd1exp(:,end),'s','DisplayName','ETD1exp');
% legend;
% hold off;


figure;
% contourf(u','LevelStep',0.1,'EdgeColor','none'); 
imagesc(uetd1exp');
colormap jet; colorbar;
set(gca,'YDir','normal');
axis square;

figure; plot(uetd1exp(1,:),'-o','LineWidth',1);

%% wavenumber
figure; hold on;
plot(k_vec,abs(fft(uetd1exp(:,end))),'-o','LineWidth',1);
xlim([min(k_vec) max(k_vec)]);
set(gca,'TickLabelInterpreter','tex','FontSize',15);
axis square;
box("on");
xlabel('$k$','Interpreter','latex','FontSize',30);
ylabel('$\hat{u}_{ss}$','Interpreter','latex','FontSize',30);

%%
figure; hold on;

npad = 10000; % # of trailing zeros u is to be padded with 
             % for fft calculation
Y = fft(u(:,end),npad);
P = abs(Y/sqrt(npad)).^2;
f = spatial_res * ((-npad/2):(npad/2-1))/npad;
f = fftshift(f);
% plot(f,P(1:npad/2+1),'-o');
plot(f,abs(Y),'-o');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
axis square;
box("on");
xlabel('$k$','Interpreter','latex','FontSize',30);
ylabel('$|\hat{u}|$','Interpreter','latex','FontSize',30);

%% figure;
figure;
npad = 10000; % # of trailing zeros u is to be padded with
% for fft calculation
f = spatial_res * ((-npad/2):(npad/2-1))/npad;
f = fftshift(f);
Yse = fft(u(:,end),npad);
Yse = abs(Yse/sqrt(npad)).^2;
plot(f,Yse,'-s','DisplayName','Spectral');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
axis square;
box("on");
xlabel('$k$','Interpreter','latex','FontSize',30);
ylabel('$|\hat{u}|$','Interpreter','latex','FontSize',30);
legend;

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