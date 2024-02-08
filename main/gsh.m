%% GSH
% This script simulates the generalized Swift Hohenberg equation (GSH).

% now using parameter array as a struct - better usability in function and
% no dependence on where a variable in the array
% -- Aditya Raj, Jan 31, 2024

close all;
clear all;

addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

%% Video switch
makeVideo = 0;

%% control parameters
eps = 0.7;

% ----mean flow----
sig = 2;
csq = 0.1; c = sqrt(csq);
gm = 50;

% ----no mean flow----
% sig = 0;
% c = 0;
% gm = 0;

% ----discretization----
lam_0 = 2*pi; % critical wavelength
dx = lam_0/8; % node spacing
dy = dx;
delta = dx;

dt = 0.2; % time step size
h = dt;

time_units = 100;
tmax = time_units/dt; % total time steps

% ----domain size----
rolls = 5;
Nx = round(rolls*lam_0/dx); % spatial nodes in x
Ny = round(rolls*lam_0/dy); % spatial nodes in y
% Nx = 5;
% Ny = 5;
N = Nx*Ny;

gmPos = 2; % position of gm; 1 for eqn 1, 2 for eqn 2 (Karimi)
gsiters = 5; % number of iterations for gauss seidel

% ----random seed----
seed = 1;
rng(seed,"twister");

% ------- running filename ---------
run_name = join([num2str(Nx) 'x' num2str(Ny) 'eps' num2str(round(100*eps)) ...
            'sig' num2str(sig) 'csq0-' num2str(round(100*csq)) ...
            'gm' num2str(gm) 't' num2str(time_units) 's' num2str(seed)]);

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


% ---- parameter struct(ure) ----
p = struct('eps', eps, 'sig', sig, 'c', c, 'gm', gm, 'dt', dt, ...
    'dx', dx, 'dy', dy, 'tmax', tmax, 'Nx', Nx, 'Ny', Ny, 'seed', seed, ...
    'I', I, 'J', J, 'Ip1', Ip1, 'Im1', Im1,'Ip2', Ip2,'Im2', Im2,'Jp1', ...
    Jp1,'Jm1', Jm1,'Jp2', Jp2,'Jm2', Jm2, 'gmPos', gmPos);

%% Initial conditions and preallocation

psimat = zeros(Nx,Ny,1); % scalar field
psivec = zeros(N,1);
% psilat(:,:,1) = -1+(1+1)*rand(Nx,Ny);
psimat(:,:,1) = -0.1 + (0.1+0.1)*rand(Nx,Ny);
psivec(:,1) = latToVec(psimat(:,:,1));


% -------- periodic boundary IC for zeta and omega ---------
% zetamat = zeros(Nx,Ny,tmax); % scalar field
% zetavec = zeros(N,tmax);
% 
% avec = 1:N;
% amat = vecToLat(avec, Nx, Ny);
% zetamat(:,:,1) = sin(amat).*cos(amat);
% zetavec(:,1) = latToVec(zetamat(:,:,1));
% 
% ulat = zeros(Nx,Ny,tmax);
% vlat = zeros(Nx,Ny,tmax);
% 
% omzmat = zeros(Nx,Ny,1); % scalar field
% omzvec = zeros(N,1);
% omzmat(:,:,1) = - laplacian(zetamat(:,:,1),para) + 0.1*rand(Nx,Ny);
% omzvec(:,1) = latToVec(omzmat(:,:,1));


%-------------- random IC for zeta and omega -----------------
zetamat = zeros(Nx,Ny,tmax); % scalar field
zetavec = zeros(N,tmax);
% omzmat = -0.01 + (0.01+0.01)*rand(Nx,Ny); % scalar field

rng(seed + 1,"twister");
% omzmat = -0.01 + (0.01+0.01) * rand(Nx,Ny);
omzmat = -1 + (1+1)*rand(Nx,Ny);
omzvec = zeros(N,1);

tic
fprintf('Constructing coefficient matrices...\n');

%% coefficients

% zeta
z1 = 1/dx^2;
z2 = -2 * ((1/dx^2) + (1/dy^2));
z3 = 1/dy^2;


a0 = - (eps-1)*h/2 + 10*h/delta^4 - 4*h/delta^2;

a1 = 1 + a0;
a2 = -4*h/delta^4 + h/delta^2;
a3 = h/delta^4;
a4 = h/(2*delta^4);

b1 = 1 - a0;
b2 = -a2;
b3 = -a3;
b4 = -a4;

% omz
m0 = sig*h/(dx^2) + sig*h/(dy^2) + h*sig*(c^2)/2;
m1 = 1 + m0;
m2 = - sig*h/(2*dx^2);
m3 = - sig*h/(2*dy^2);
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

a = zeros(Nx,Ny);
b = zeros(Nx,Ny);
A = zeros(N,N);
B = zeros(N,N);

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

        a(I(i),J(j)) = a1;

        a(I(i),Jp1(j)) = a2;
        a(I(i),Jm1(j)) = a2;
        a(Ip1(i),J(j)) = a2;
        a(Im1(i),J(j)) = a2;

        a(Ip1(i),Jp1(j)) = a3;
        a(Ip1(i),Jm1(j)) = a3;
        a(Im1(i),Jp1(j)) = a3;
        a(Im1(i),Jm1(j)) = a3;

        a(Im2(i),J(j)) = a4;
        a(I(i),Jp2(j)) = a4;
        a(I(i),Jm2(j)) = a4;
        a(Ip2(i),J(j)) = a4;

        A(n,:) = latToVec(a)';

        b = -1*a;

        b(I(i),J(j)) = b1;

        B(n,:) = latToVec(b)';
        
        ml(I(i),J(j)) = m1;
        ml(Im1(i),J(j)) = m2;
        ml(Ip1(i),J(j)) = m2;
        ml(I(i),Jm1(j)) = m3;
        ml(I(i),Jp1(j)) = m3;

        Ml(n,:) = latToVec(ml)';

        mr = - ml;

        mr(I(i),J(j)) = 1 - m0;

        Mr(n,:) = latToVec(mr)';

        n = n+1;
        z = zeros(Nx,Ny);
        a = zeros(Nx,Ny);
        b = zeros(Nx,Ny);
        ml = zeros(Nx,Ny);
        mr = zeros(Nx,Ny);
    end
end

% Z = sparse(Z);

Lz = tril(Z);
Uz = triu(Z,1);
Lzinv = inv(Lz);

matdivpsi = A\B;
matdivomz = Ml\Mr;
zetavect = rand(N,1);
Zzeta = zeros(N,1);

 toc

%% Time integration
fprintf('Running time integration...\n');

for t = 1:tmax    
    %------------------ zeta, iterative ------------------------
    zetamat(:,:,t) = zetafunc(omzmat(:,:,t),zetamat(:,:,t),p);
    [ulat(:,:,t),vlat(:,:,t)] = uvzeta(zetamat(:,:,t),p);
    
    %------------------ explicit nonlinear --------------------
    psitilde = rk2eq1(psimat(:,:,t),zetamat(:,:,t),p,@nlpgsh1); 
    psitilde = latToVec(psitilde);

    omztilde = rk2eq2(omzmat(:,:,t),psimat(:,:,t),p,@nlpgsh2);
    omztilde = latToVec(omztilde);
    
    %------------------ implicit CN linear ----------------------
    psivec(:,t+1) = matdivpsi * psitilde;
    omzvec(:,t+1) = matdivomz * omztilde;

    %--------- converting from vector to matrix ----------
    psimat(:,:,t+1) = vecToLat(psivec(:,t+1),Nx,Ny);
    omzmat(:,:,t+1) = vecToLat(omzvec(:,t+1),Nx,Ny);
    
    zetamat(:,:,t+1) = zetamat(:,:,t);
    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
end
toc

%% dynamics video

if makeVideo == 1
    dynvideoname = 'gshDyn';
    dynVideoFilename = join([dynvideoname run_name]);
    lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
    %lat_dyn_video.FrameRate = 30;
    open(lat_dyn_video);
    [X,Y] = meshgrid(1:Nx,1:Ny);
    figure();
    % hold on;
    for t = 1:altframes:tmax
        % imagesc(psilat(:,:,t));
        plot1 = contourf(psimat(:,:,t),'levels',0.1, 'Linecolor', 'none');
        set(gca,'YDir','normal');
        hold on;
        % imagesc(zetalat(:,:,t));
        colorbar;
        colormap jet
        clim([-0.8 0.8]);
        plot2 = quiver(X,Y,vlat(:,:,t),ulat(:,:,t),2,'black');
        xlim([1 Nx]);
        ylim([1 Ny]);
        frame = getframe(gcf);
        writeVideo(lat_dyn_video,frame);
        clear plot1 plot2;
        hold off;
    end
    hold off;
    close(lat_dyn_video);
end

%% psi figure

% figure;
% imagesc(psimat(:,:,tmax)');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;
% title('\psi');

% figure;
% imagesc(ulat(:,:,6)');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;

%%
figure;
[X,Y] = meshgrid(1:Nx,1:Ny);
contourf(psimat(:,:,tmax),'levels',0.1, 'Linecolor', 'none');
clim([-1 1]);
hold on;
quiver(X,Y,vlat(:,:,tmax),ulat(:,:,tmax),2,'black');
hold off;
set(gca,'YDir','normal');
xlim([1 Nx]);
ylim([4 Ny]);
colorbar;
colormap jet;
title('\psi');

%% zeta contour figure

figure;
contourf(zetamat(:,:,end),'Linecolor','none');
set(gca,'YDir','normal');
colormap jet;
colorbar;
title('\zeta');


%% zeta contour figure with rolls
timestep = tmax; % change this to desired time step (where spirals are seen)
figure;
p1 = contourf(abs(zetamat(:,:,timestep)),'levels',0.1, 'Linecolor', 'none');
set(gca,'YDir','normal');
xlim([1 Nx]);
ylim([1 Ny]);
colormap jet;
colorbar;
hold on;
p2 = quiver(X,Y,2*vlat(:,:,timestep),2*ulat(:,:,timestep),2,'black');
% clim([0 5]);
p3 = contourf(psimat(:,:,timestep), 'levels', 1, 'Linecolor', 'black', ...
    'Linewidth', 3,'Facecolor', 'none');
hold off;

%% Functions

function zeta = zetafunc(omzmat, zetamat, p)
    
    delta = p.dx;

    I = p.I;
    J = p.J;
    Ip1 = p.Ip1;
    Im1 = p.Im1;
    Jp1 = p.Jp1;
    Jm1 = p.Jm1;

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


function lapl = laplacian(var, p)

    delta = p.dx;

    I = p.I;
    J = p.J;
    Ip1 = p.Ip1;
    Im1 = p.Im1;
    Jp1 = p.Jp1;
    Jm1 = p.Jm1;
    
    lapl = zeros(size(var));

    for i = 1:length(I)
        for j = 1:length(J)
            lapl(i,j) = (1/(delta^2)) * ...
                        ( var(Im1(i),J(j)) + var(Ip1(i),J(j)) ...
                            - 4 * var(I(i),J(j)) + ...
                        + var(I(i),Jm1(j)) + var(I(i),Jp1(j)) );
        end
    end
end

function [u,v] = uvzeta(zeta,p)
    
    dx = p.dx;
    dy = p.dy;
    Nx = p.Nx;
    Ny = p.Ny;

    I = p.I;
    J = p.J;
    Ip1 = p.Ip1;
    Im1 = p.Im1;
    Jp1 = p.Jp1;
    Jm1 = p.Jm1;

    u = zeros(Nx,Ny);
    v = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            dzetadx = (zeta(Ip1(i),J(j)) - zeta(Im1(i),J(j)))/(2*dx);
            dzetady = (zeta(I(i),Jp1(j)) - zeta(I(i),Jm1(j)))/(2*dy);
            u(i,j) = dzetady;
            v(i,j) = -dzetadx;
        end
    end
end

function np1 = nlpgsh1(psit,zetat,p)    
    dx = p.dx;
    dy = p.dy;
    Nx = p.Nx;
    Ny = p.Ny;
    gmPos = p.gmPos;

    I = p.I;
    J = p.J;
    Ip1 = p.Ip1;
    Im1 = p.Im1;
    Jp1 = p.Jp1;
    Jm1 = p.Jm1;

    % c = 1/(4*dx*dy);
    np1 = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            dzetadx = (zetat(Ip1(i),J(j)) - zetat(Im1(i),J(j)))/(2*dx);
            dzetady = (zetat(I(i),Jp1(j)) - zetat(I(i),Jm1(j)))/(2*dy);
            dpsidx = (psit(Ip1(i),J(j)) - psit(Im1(i),J(j)))/(2*dx);
            dpsidy = (psit(I(i),Jp1(j)) - psit(I(i),Jm1(j)))/(2*dy);
            psi = psit(I(i),J(j));
            % if gmPos == 1
            %     np1(i,j) = - 3*psi^3 - gm* ( dzetady * dpsidx - dzetadx * dpsidy );
            % elseif gmPos == 2
                np1(i,j) = - psi^3 - dzetady * dpsidx + dzetadx * dpsidy ;
            % else
            %     error('gmPos must be either 1 or 2');
            % end
        end
    end
end


function np2 = nlpgsh2(omz,psi,p)   
    gm = p.gm;
    dx = p.dx;
    dy = p.dy;
    Nx = p.Nx;
    Ny = p.Ny;
    gmPos = p.gmPos;

    I = p.I;
    J = p.J;
    Ip1 = p.Ip1;
    Im1 = p.Im1;
    Jp1 = p.Jp1;
    Jm1 = p.Jm1;
    
    Ip2 = p.Ip2;
    Im2 = p.Im2;
    Jp2 = p.Jp2;
    Jm2 = p.Jm2;

    np2 = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            dpsidx = psi(Ip1(i),J(j))/(2*dx) - psi(Im1(i),J(j))/(2*dx);
            dpsidy = psi(I(i),Jp1(j))/(2*dy) - psi(I(i),Jm1(j))/(2*dy);
            d3psidx3 = (1/(2*dx^3)) * (-psi(Im2(i),J(j)) + 2 * psi(Im1(i),J(j)) ...
                    - 2 * psi(Ip1(i),J(j)) + psi(Ip2(i),J(j)));
            d3psidy3 = (1/(2*dy^3)) * (-psi(I(i),Jm2(j)) + 2 * psi(I(i),Jm1(j)) ...
                    - 2 * psi(I(i),Jp1(j)) + psi(I(i),Jp2(j)));
            d3psidxdy2 = (1/(2*dx*dy^2)) * (psi(Ip1(i),Jm1(j)) - psi(Im1(i),Jm1(j)) ...
                    - 2 * psi(Ip1(i),J(j)) + 2*psi(Im1(i),J(j)) ...
                    + psi(Ip1(i),Jp1(j)) - psi(Im1(i),Jp1(j)));
            d3psidydx2 = (1/(2*dy*dx^2)) * (psi(Im1(i),Jp1(j)) - psi(Im1(i),Jm1(j)) ...
                    - 2 * psi(I(i),Jp1(j)) + 2 * psi(I(i),Jm1(j)) ...
                    + psi(Ip1(i),Jp1(j)) - psi(Ip1(i),Jm1(j)));
            % if gmPos == 1
            %     np2(i,j) = ( dpsidy * ( d3psidx3 + d3psidxdy2 ) ...
            %         - dpsidx * ( d3psidydx2 + d3psidy3 ) );
            % elseif gmPos == 2
            np2(i,j) = -gm * ( dpsidy * (d3psidx3 + d3psidxdy2) ...
                    - dpsidx * (d3psidydx2 + d3psidy3) );
            % else
            %     error('gmPos must be either 1 or 2');
            % end
        end
    end
end


function psitilde = rk2eq1(psi,zeta,p,dynfunc)
    dt = p.dt;
    k1 = dynfunc(psi,zeta,p);
    %u1 = u+k1*dt;
    k2 = dynfunc(psi+k1*dt,zeta,p);
    psitilde = psi + dt*((k1+k2)/2);
end

function omztilde = rk2eq2(omz,psi,p,dynfunc)
    dt = p.dt;
    k1 = dynfunc(omz,psi,p);
    %u1 = u+k1*dt;
    k2 = dynfunc(omz+k1*dt,psi,p);
    omztilde = omz + dt*((k1+k2)/2);
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


function xnew = gauss_seidel(A, b, x)
    N = length(x);
    xnew = zeros(N,1);
    iters = 0;
    while max(abs(A*x-b)) > 1e-5
        for i = 1:N
            xnew(i) = ( b(i) - ( sum(A(i,:)*x) - A(i,i)*x(i) ) )/A(i,i);
            x(i) = xnew(i);
        end
        iters = iters+1;
    end
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
