function p = paramgsh(eps,sig,csq,gm,dx,dy,dt,rolls,trtimeu,totimeu,seed)
% p = paramgsh(eps,sig,csq,gm,dx,dy,dt,rolls,timeUnits,seed)
% lam_0 is 2*pi
% dx or dy should be taken as 2*pi/k, k is an integer


%----epsilon------
% eps = 0.7;

% ----mean flow----
% sig = 2;
% csq = 0.1; 
% 
c = sqrt(csq);
% gm = 50;

% ----no mean flow----
% sig = 0;
% c = 0;
% gm = 0;

% ----discretization----
lam_0 = 2*pi; % critical wavelength
% dx = lam_0/8; % node spacing
% dy = dx;
% delta = dx;

% dt = 0.1; % time step size
% h = dt;

% timeUnits = 25;
tmax = totimeu/dt; % total time steps


% ----domain size----
% rolls = 5;
Nx = round(rolls*lam_0/dx); % spatial nodes in x
Ny = round(rolls*lam_0/dy); % spatial nodes in y
% Nx = 5;
% Ny = 5;
N = Nx*Ny;

gmPos = 2; % position of gm; 1 for eqn 1, 2 for eqn 2 (Karimi)
gsiters = 5; % number of iterations for gauss seidel

% ----random seed----
% seed = 1;
rng(seed,"twister");

% ------- running filename ---------
run_name = join([num2str(Nx) 'x' num2str(Ny) 'eps' num2str(round(100*eps)) ...
            'sig' num2str(sig) 'csq0-' num2str(round(100*csq)) ...
            'gm' num2str(gm) 't' num2str(totimeu) 's' num2str(seed)]);

% ----video----
altframes = 2;

% ------ vectors for periodic BCs -----------
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


% ----  (p)arameter (struct)ure ----
p = struct('eps', eps, 'sig', sig, 'c', c, 'gm', gm, 'dt', dt, ...
    'dx', dx, 'dy', dy, 'tmax', tmax, 'Nx', Nx, 'Ny', Ny, 'N',N, ...
    'seed', seed, 'I', I, 'J', J, 'Ip1', Ip1, 'Im1', Im1,'Ip2', Ip2, ...
    'Im2', Im2,'Jp1', Jp1,'Jm1', Jm1,'Jp2', Jp2,'Jm2', Jm2, ...
    'gmPos', gmPos, 'run_name', run_name, 'lam_0', lam_0, ...
    'trtimeu', trtimeu, 'totimeu', totimeu, 'rolls', rolls);

end

