%% GSH with TS calculations
% This script simulates the generalized Swift Hohenberg equation (GSH) and
% leading Lyapunov exponents. (in progress)
% -- Aditya Raj, Jan 31, 2024

close all;
clear all;

addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

%% switches
makeDynVideo = 0;

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

dt = 0.05; % time step size
h = dt;

time_units = 25;
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


% ----  (p)arameter (struct)ure ----
p = struct('eps', eps, 'sig', sig, 'c', c, 'gm', gm, 'dt', dt, ...
    'dx', dx, 'dy', dy, 'tmax', tmax, 'Nx', Nx, 'Ny', Ny, 'seed', seed, ...
    'I', I, 'J', J, 'Ip1', Ip1, 'Im1', Im1,'Ip2', Ip2,'Im2', Im2,'Jp1', ...
    Jp1,'Jm1', Jm1,'Jp2', Jp2,'Jm2', Jm2, 'gmPos', gmPos);

% p = struct2table(p);

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
seed = 4;
rng(seed,"twister");
% omzmat = -0.01 + (0.01+0.01) * rand(Nx,Ny);
omzmat = -1 + (1+1)*rand(Nx,Ny);
omzvec = zeros(N,1);

%-------------- tangent space psi, zeta and omega -----------------
dpsi1 = zeros(Nx,Ny,tmax);
dpsi1vec = zeros(N,1);
domz1 = zeros(Nx,Ny,tmax);
domz1vec = zeros(N,1);
dzeta1 = zeros(Nx,Ny,tmax);

dpsi1(:,:,1) = rand(Nx,Ny);
dpsi1(:,:,1) = dpsi1(:,:,1)./norm(latToVec(dpsi1(:,:,1)));

domz1(:,:,1) = rand(Nx,Ny);
domz1(:,:,1) = domz1(:,:,1)./norm(latToVec(domz1(:,:,1)));

fpv = zeros(2*N,1);
fpv(:,1) = [latToVec(dpsi1(:,:,1)); latToVec(domz1(:,:,1))];
fpv(:,1) = [latToVec(dpsi1(:,:,1)); latToVec(domz1(:,:,1))]./ norm(fpv(:,1));
fpvmag = norm(fpv(:,1));

%% lambda1 calculation preallocation
normtime = 2;
lam1inst = [];



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

% TS 2
e0 = -m0;
e1 = 1+e0;
e2 = -m2;
e3 = -m3;


%% coefficient matrix construction

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
    psitilde = rk2gsh1(psimat(:,:,t),zetamat(:,:,t),p,@nlpgsh1); 
    psitilde = latToVec(psitilde);

    omztilde = rk2gsh2(omzmat(:,:,t),psimat(:,:,t),p,@nlpgsh2);
    omztilde = latToVec(omztilde);
    
    %------------------ implicit CN linear ----------------------
    psivec(:,t+1) = matdivpsi * psitilde;
    omzvec(:,t+1) = matdivomz * omztilde;

    %--------- converting from vector to matrix ----------
    psimat(:,:,t+1) = vecToLat(psivec(:,t+1),Nx,Ny);
    omzmat(:,:,t+1) = vecToLat(omzvec(:,t+1),Nx,Ny);
    
    %------------------ explicit TS linear ----------------------
    dzeta1(:,:,t) = zetafunc(domz1(:,:,t),dzeta1(:,:,t),p);
    dpsi1(:,:,t+1) = rk2ts1(psimat(:,:,t),zetamat(:,:,t),dpsi1(:,:,t), ...
                        dzeta1(:,:,t),p,@ts1exp,@fds);
    domz1(:,:,t+1) = rk2ts2(omzmat(:,:,t),psimat(:,:,t),domz1(:,:,t), ...
                        dpsi1(:,:,t),p,@ts2exp,@fds);
    fpv(:,t+1) = [latToVec(dpsi1(:,:,t)); latToVec(domz1(:,:,t))];
    fpvmag(t) = norm(fpv(:,t));
    
    %------------------ renormalization and LLE ----------------------
    if rem(t,normtime) == 0
        lam1inst = [lam1inst, log(abs(norm(fpv(:,t+1))/norm(fpv(:,t))))];
        fpv(:,t+1) = fpv(:,t+1)./norm(fpv(:,t+1));
        dpsi1(:,:,t+1) = vecToLat(fpv(1:N,t+1),Nx,Ny);
        domz1(:,:,t+1) = vecToLat(fpv(N+1:2*N,t+1),Nx,Ny);
    end
    
    zetamat(:,:,t+1) = zetamat(:,:,t);
    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
    
end
toc


%% dynamics video

if makeDynVideo == 1
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
contourf(psimat(:,:,tmax),'levels',0.1, 'Linecolor', 'none');
clim([-1 1]);
hold on;
[X,Y] = meshgrid(1:Nx,1:Ny);
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
% timestep = tmax; % change this to desired time step (where spirals are seen)
% figure;
% p1 = contourf(abs(zetamat(:,:,timestep)),'levels',0.1, 'Linecolor', 'none');
% set(gca,'YDir','normal');
% xlim([1 Nx]);
% ylim([1 Ny]);
% colormap jet;
% colorbar;
% hold on;
% p2 = quiver(X,Y,2*vlat(:,:,timestep),2*ulat(:,:,timestep),2,'black');
% % clim([0 5]);
% p3 = contourf(psimat(:,:,timestep), 'levels', 1, 'Linecolor', 'black', ...
%     'Linewidth', 3,'Facecolor', 'none');
% hold off;

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

function dpsirhs = ts1exp(psi,zeta,dpsi,dzeta,p,fdsfunc)
    Nx = p.Nx;
    Ny = p.Ny;

    eps = p.eps;

    I = p.I;
    J = p.J;

    dpsirhs = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            dzetady = fdsfunc(zeta,0,1,i,j,p);
            dzetadx = fdsfunc(zeta,1,0,i,j,p);
            dpsidy = fdsfunc(psi,0,1,i,j,p);
            dpsidx = fdsfunc(psi,1,0,i,j,p);

            ddzetady = fdsfunc(dzeta,0,1,i,j,p);
            ddzetadx = fdsfunc(dzeta,1,0,i,j,p);
            ddpsidy = fdsfunc(dpsi,0,1,i,j,p);
            ddpsidx = fdsfunc(dpsi,1,0,i,j,p);
            
            d4dpsidx4 = fdsfunc(dpsi,4,0,i,j,p);
            d4dpsidx2dy2 = fdsfunc(dpsi,2,2,i,j,p);
            d4dpsidy4 = fdsfunc(dpsi,0,4,i,j,p);

            d2dpsidx2 = fdsfunc(dpsi,2,0,i,j,p);
            d2dpsidy2 = fdsfunc(dpsi,0,2,i,j,p);

            dpsirhs(i,j) = - dzetady * ddpsidx + dzetadx * ddpsidy ...
                       - ddzetady * dpsidx + ddzetadx * dpsidy ...
                       - d4dpsidx4 - 2 * d4dpsidx2dy2 - d4dpsidy4 ...
                       - 2 * d2dpsidx2 - 2 * d2dpsidy2 ...
                       + (eps - 1 - 3 * psi(I(i),J(j))^2) * dpsi(I(i),J(j));
        end
    end
end

function domzrhs = ts2exp(omz,psi,domz,dpsi,p,fdsfunc)
    Nx = p.Nx;
    Ny = p.Ny;
    
    sig = p.sig;
    c = p.c;
    gm = p.gm;

    I = p.I;
    J = p.J;

    domzrhs = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            
            d2domzdx2 = fdsfunc(domz,2,0,i,j,p);
            d2domzdy2 = fdsfunc(domz,0,2,i,j,p);

            dpsidy = fdsfunc(psi,0,1,i,j,p);
            dpsidx = fdsfunc(psi,1,0,i,j,p);

            ddpsidy = fdsfunc(dpsi,0,1,i,j,p);
            ddpsidx = fdsfunc(dpsi,1,0,i,j,p);
            
            d3dpsidx3 = fdsfunc(dpsi,3,0,i,j,p);
            d3dpsidxdy2 = fdsfunc(dpsi,1,2,i,j,p);
            d3dpsidydx2 = fdsfunc(dpsi,2,1,i,j,p);
            d3dpsidy3 = fdsfunc(dpsi,0,3,i,j,p);
            
            d3psidx3 = fdsfunc(psi,3,0,i,j,p);
            d3psidxdy2 = fdsfunc(psi,1,2,i,j,p);
            d3psidydx2 = fdsfunc(psi,2,1,i,j,p);
            d3psidy3 = fdsfunc(psi,0,3,i,j,p);

            domzrhs(i,j) = sig * d2domzdx2 + sig * d2domzdy2 ...
                        - (c^2) * domz(I(i),J(j)) ...
                       - gm * dpsidy * d3dpsidx3 - gm * dpsidy * d3dpsidxdy2 ...
                       + gm * dpsidx * d3dpsidydx2 + gm * dpsidx * d3dpsidy3 ...
                       - gm * ddpsidy * d3psidx3 - gm * ddpsidy * d3psidxdy2 ...
                       + gm * ddpsidx * d3psidydx2 + gm * ddpsidx * d3psidy3 ;
        end
    end
end


function psitilde = rk2gsh1(psi,zeta,p,dynfunc)
    dt = p.dt;
    k1 = dynfunc(psi,zeta,p);
    %u1 = u+k1*dt;
    k2 = dynfunc(psi+k1*dt,zeta,p);
    psitilde = psi + dt*((k1+k2)/2);
end

function omztilde = rk2gsh2(omz,psi,p,dynfunc)
    dt = p.dt;
    k1 = dynfunc(omz,psi,p);
    %u1 = u+k1*dt;
    k2 = dynfunc(omz+k1*dt,psi,p);
    omztilde = omz + dt*((k1+k2)/2);
end

function dpsinp1 = rk2ts1(psi,zeta,dpsi,dzeta,p,tsdynfunc,fdsfunc)
    dt = p.dt;
    k1 = tsdynfunc(psi,zeta,dpsi,dzeta,p,fdsfunc);
    %u1 = u+k1*dt;
    k2 = tsdynfunc(psi,zeta,dpsi+k1*dt,dzeta,p,fdsfunc);
    dpsinp1 = dpsi + dt*((k1+k2)/2);
end

function domznp2 = rk2ts2(omz,psi,domz,dpsi,p,tsdynfunc,fdsfunc)
    dt = p.dt;
    k1 = tsdynfunc(omz,psi,domz,dpsi,p,fdsfunc);
    %u1 = u+k1*dt;
    k2 = tsdynfunc(omz,psi,domz+k1*dt,dpsi,p,fdsfunc);
    domznp2 = domz + dt*((k1+k2)/2);
end



function derivtv = fds(v,xd,yd,i,j,p) 
% dvxddyd = fds(v,xd,yd,i,j,p) 
% Yields the spatial finite difference derivative
% Not all derivatives are possible.
% Example: d3v/dxdy2(i,j) = fds(v,1,2,i,j,p)
    dx = p.dx;
    dy = p.dy;

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

    % for i = 1:length(I)
    %     for j = 1:length(J)

    if xd == 1 && yd == 0 % dv/dx
        derivtv = v(Ip1(i),J(j))/(2*dx) - v(Im1(i),J(j))/(2*dx);
    elseif xd == 0 && yd == 1 % dv/dy
        derivtv = v(I(i),Jp1(j))/(2*dy) - v(I(i),Jm1(j))/(2*dy);
    elseif xd == 2 && yd == 0
        derivtv = (1/(dx^2)) * (v(Im1(i),J(j)) - 2 * v(I(i),J(j)) ...
                                    + v(Ip1(i),J(j)));
    elseif xd == 0 && yd == 2
        derivtv = (1/(dy^2)) * (v(I(i),Jm1(j)) - 2 * v(I(i),J(j)) ...
                                    + v(I(i),Jp1(j)));
    elseif xd == 3 && yd == 0 % d3v/dx3
        derivtv = (1/(2*dx^3)) * (-v(Im2(i),J(j)) + 2 * v(Im1(i),J(j)) ...
            - 2 * v(Ip1(i),J(j)) + v(Ip2(i),J(j)));
    elseif xd == 0 && yd == 3 % d3v/dy3
        derivtv = (1/(2*dy^3)) * (-v(I(i),Jm2(j)) + 2 * v(I(i),Jm1(j)) ...
            - 2 * v(I(i),Jp1(j)) + v(I(i),Jp2(j)));
    elseif xd == 1 && yd == 2 % d3v/dxdy2
        derivtv = (1/(2*dx*dy^2)) * (v(Ip1(i),Jm1(j)) - v(Im1(i),Jm1(j)) ...
            - 2 * v(Ip1(i),J(j)) + 2*v(Im1(i),J(j)) ...
            + v(Ip1(i),Jp1(j)) - v(Im1(i),Jp1(j)));
    elseif xd == 2 && yd == 1 % d3v/dydx2
        derivtv = (1/(2*dy*dx^2)) * (v(Im1(i),Jp1(j)) - v(Im1(i),Jm1(j)) ...
            - 2 * v(I(i),Jp1(j)) + 2 * v(I(i),Jm1(j)) ...
            + v(Ip1(i),Jp1(j)) - v(Ip1(i),Jm1(j)));
    elseif xd == 4 && yd == 0 %d4v/dx4
        derivtv = (1/(dx^4)) * (v(Im2(i),J(j)) - 4 * v(Im1(i),J(j)) ...
            + 6*v(I(i),J(j)) - 4 * v(Ip1(i),J(j)) + v(Ip2(i),J(j)));
    elseif xd == 0 && yd == 4 % d4v/dy4
        derivtv = (1/(dy^4)) * (v(I(i),Jm2(j)) - 4 * v(I(i),Jm1(j)) ...
            + 6*v(I(i),J(j)) - 4 * v(I(i),Jp1(j)) + v(I(i),Jp2(j)));
    elseif xd == 2 && yd == 2 % d4v/dx2dy2
        derivtv = (1/(dx^2*dy^2)) * ( v(Ip1(i),Jp1(j)) + v(Ip1(i),Jm1(j))  ...
                    + v(Im1(i),Jp1(j)) + v(Im1(i),Jm1(j))) ...
                - (2/(dx^2*dy^2)) * ( v(Ip1(i),J(j)) + v(Im1(i),J(j)) ...
                    + v(I(i),Jm1(j)) +  v(I(i),Jp1(j)) ) ...
                + (4/(dx^2*dy^2)) * ( v(I(i),J(j)) );
    else
        error(['Requested spatial finite difference ' ...
            'derivative is not available!']);
    end
    %     end
    % end
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


% function u = rk4(u,dynFunc,h)
%     % k values
%     k1 = dynFunc(u);
%     k2 = dynFunc(u + (0.5*h)*k1);
%     k3 = dynFunc(u + (0.5*h)*k2);
%     k4 = dynFunc(u + h*k3);
%     % dynamics
%     u = u + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
% end
% 
% 
% function xnew = gauss_seidel(A, b, x)
%     N = length(x);
%     xnew = zeros(N,1);
%     iters = 0;
%     while max(abs(A*x-b)) > 1e-5
%         for i = 1:N
%             xnew(i) = ( b(i) - ( sum(A(i,:)*x) - A(i,i)*x(i) ) )/A(i,i);
%             x(i) = xnew(i);
%         end
%         iters = iters+1;
%     end
% end


% function mat = jacob(psi,omz,zeta,t,p,fdsfunc)
% % under construction, do not use
%     dx = p.dx;
%     dy = p.dy;
% 
%     I = p.I;
%     J = p.J;
%     Ip1 = p.Ip1;
%     Im1 = p.Im1;
%     Jp1 = p.Jp1;
%     Jm1 = p.Jm1;
% 
%     Ip2 = p.Ip2;
%     Im2 = p.Im2;
%     Jp2 = p.Jp2;
%     Jm2 = p.Jm2;
% 
%     a = zeros(Nx,Ny);
%     b = zeros(Nx,Ny);
% 
%     m0 = -3*dt/dx^4 + 4*dt/(dx^2*dy^2) - 3*dt/dy^4 + 2*dt/dx^2 + 2*dt/dy^2;
%     a2b = 2*dt/dx^4 + 2*dt/(dx^2*dy^2) - dt/dx^2;
%     a4b = 2*dt/dy^4 + 2*dt/(dx^2*dy^2) - dt/dy^2;
% 
%     a6 = -dt/(dx^2*dy^2);
%     a7 = -dt/(2*dx^4);
%     a8 = -dt/(2*dy^4);
% 
%     for i = 1:length(I)
%         for j = 1:length(J)
%             a0 = m0 + (dt/2)*(eps-1-3*psi(I(i),J(j))^2);
%             a1 = 1+a0;
%             a2a = -(dt/(4*dx)) * fdsfunc(zeta,0,1,p);
%             a2 = a2a + a2b;
%             a3 = -a2a + a2b;
%             a4a = -(dt/(4*dy)) * fdsfunc(zeta,1,0,p);
%             a4 = a4a + a4b;
%             a5 = -a4a + a4b;
%             a6 = -dt/(dx^2*dy^2);
%             a9 = (dt/(4*dx)) * fdsfunc(psi,0,1,p);
%             a10 = -a9;
%             a11 = -(dt/(4*dy)) * fdsfunc(psi,1,0,p);
%             a12 = -a11;
% 
%             a(I(i),J(j)) = a1;
% 
%             a(Ip1(i),J(j)) = a2;
%             a(Im1(i),J(j)) = a3;
%             a(I(i),Jp1(j)) = a4;
%             a(I(i),Jm1(j)) = a5;
% 
%             a(Ip1(i),Jp1(j)) = a6;
%             a(Ip1(i),Jm1(j)) = a6;
%             a(Im1(i),Jp1(j)) = a6;
%             a(Im1(i),Jm1(j)) = a6;
% 
%             a(Im2(i),J(j)) = a7;
%             a(Ip2(i),J(j)) = a7;
%             a(I(i),Jp2(j)) = a8;
%             a(I(i),Jm2(j)) = a8;
% 
%             a(I(i),J(j)) = 1 - a0;
% 
%             b(Ip1(i),J(j)) = -a2;
%             b(Im1(i),J(j)) = -a3;
%             b(I(i),Jp1(j)) = -a4;
%             b(I(i),Jm1(j)) = -a5;
% 
%             b(Ip1(i),Jp1(j)) = -a6;
%             b(Ip1(i),Jm1(j)) = -a6;
%             b(Im1(i),Jp1(j)) = -a6;
%             b(Im1(i),Jm1(j)) = -a6;
% 
%             b(Im2(i),J(j)) = -a7;
%             b(Ip2(i),J(j)) = -a7;
%             b(I(i),Jp2(j)) = -a8;
%             b(I(i),Jm2(j)) = -a8;
% 
%         end
%     end
% end
