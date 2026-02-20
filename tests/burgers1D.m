%% 1D Burgers Equation: Finite Difference + RK2
% u_t + u u_x = nu u_xx
% Periodic boundary conditions

clear; 
% close all;

%% Default Parameters
L  = 2*pi;        % domain length
Nx = 256;         % number of grid points
dx = L / Nx;

nu = 0.01;        % viscosity
alpha_B = 0.5;

dt = 1e-3;        % time step
T  = 50;           % final time
Nt = round(T/dt);

%% Parameters changed
% L  = 2*pi;        % domain length
% Nx = 32;         % number of grid points
% dx = L / Nx;
% 
% nu = 0.03;        % viscosity
% 
% dt = 1e-2;        % time step
% T  = 5;           % final time
% Nt = round(T/dt);

%% Time vector
timevec = linspace(0,T,Nt);

%% Spatial grid
x = (0:Nx-1).' * dx;

%% Initial condition (sine-based)
% u0 = 0.1*sin(4*x)+1;                      % simple sine IC
% u0 = abs(cos(0.5*x+pi/2));              % cosine wave (positive)
% u0 = sin(x) + 0.5*sin(2*x);             % richer IC
% rng(1); u0 = -1 + 2*rand(size(x));      % random initial conditions mean 0
% u0 = 2*ones(size(x));
% u0 = sin(4*x);                          % sine waves, mean 0  
% rng(1); u0 = -0.5 + 1*rand(size(x));    % random initial conditions, mean 0
rng(1); u0 = 1 + 1*rand(size(x));   % random initial conditions, mean 5.5
% rng(1); u0 = -0.5 + 0.1*rand(size(x));  % random initial conditions, mean -0.5
u = u0;


%% Preallocate
u_new = zeros(Nx,1);

ustore = [];

ustore = zeros(Nx, Nt);  
seed = 1;
%% Time integration (RK2)
for n = 1:Nt

    % Stage 1
    k1 = burgers_rhs(u, dx, nu, alpha_B);

    % Stage 2
    u_tilde = u + dt * k1;
    k2 = burgers_rhs(u_tilde, dx, nu, alpha_B);

    % RK2 update
    u = u + 0.5 * dt * (k1 + k2);

    % add noise at every time unit
    stepsperunit = round(0.05/dt);

    if mod(n,stepsperunit) == 0
        % Add noise to the solution
        rng(seed); seed = seed +1;
        u = u + 0.1*(-0.5 + 1 * rand(size(u)));
        % u = (-0.5 + 1 * rand(size(u)));
    end

    % % add noise at every time step (all goes to zero)
    % rng(seed); seed = seed +1;
    % u = u + (-0.005 + 001 * rand(size(u)));


    ustore(:,n) = u;

    %Simple visualization
    % if mod(n,200) == 0
    %     plot(x, u, 'LineWidth', 1.5);
    %     ylim([min(u0) max(u0)]);
    %     xlabel('x'); ylabel('u');
    %     title(sprintf('t = %.2f', n*dt));
    %     drawnow;
    % end
end

%%
figure; imagesc(x,timevec,ustore'); colorbar;
set(gca,'YDir','normal','FontSize',15);
xlabel('$x$','Interpreter','latex','FontSize',30);
ylabel('$t$','Rotation',0,'Interpreter','latex','FontSize',30);
% clim([-1 1]); 
colormap jet;
% title('u'); 
axis square;

%%
figure; imagesc(x,timevec,ustore); colorbar;
set(gca,'YDir','normal','FontSize',15);
xlabel('$t$','Interpreter','latex','FontSize',30);
ylabel('$x$','Rotation',0,'Interpreter','latex','FontSize',30);
xlim([0 5]);
% clim([-1 1]); 
colormap jet;
% title('u'); 
axis square;


%% --- RHS FUNCTION ---
function rhs = burgers_rhs(u, dx, nu, alpha_B)
    % Periodic finite differences

    % Periodic shifts
    u_ip = circshift(u, -1);   % u_{i+1}
    u_im = circshift(u,  1);   % u_{i-1}

    % First derivative (central)
    ux = alpha_B * (u_ip - u_im) / (2*dx);

    % Second derivative (central)
    uxx = (u_ip - 2*u + u_im) / dx^2;

    % Burgers RHS
    rhs = -u .* ux + nu * uxx;
end
