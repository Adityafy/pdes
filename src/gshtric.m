function [psimat, omzmat, zetamat, ulat, vlat] = gshtric(p)
% [psimat, omzmat, zetamat] = gshtric(p)
% The time integration function that yields psi, omega_z, zeta as matrices 
% evolved from initial states through time amount of transient time given 
% (trtimeu). These matrices serve as the initial conditions for time
% evolution (totimeu evolution) of psi, omz, zeta that is suited for
% tangent space calculations.
addpath('../src/');

dt = p.dt; % time step size

trtimeu = p.trtimeu;
tmax = trtimeu/dt; % total time steps

Nx = p.Nx;
Ny = p.Ny;
N = Nx*Ny;

seed = p.seed;
rng(seed,"twister");


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
zetamat = zeros(Nx,Ny,1); % scalar field
ulat = zeros(Nx,Ny,1);
vlat = zeros(Nx,Ny,1);
% zetavec = zeros(N,tmax);
% omzmat = -0.01 + (0.01+0.01)*rand(Nx,Ny); % scalar field

rng(seed + 1,"twister");
% omzmat = -0.01 + (0.01+0.01) * rand(Nx,Ny);
omzmat = -1 + (1+1)*rand(Nx,Ny);
omzvec = zeros(N,1);


%% Time integration
run_name = p.run_name;
fprintf(join(['Running transients ',run_name]));

fprintf('\nCalculating coefficient matrices for implicit calcs...\n');
[matdivpsi, matdivomz] = impMatGSH(p);
toc

fprintf('Running time integration...\n');
for t = 1:tmax    
    %------------------ zeta, iterative ------------------------
    zetamat = iterativeZetaOmz(omzmat,zetamat,p);
    [ulat,vlat] = uvzeta(zetamat,p);
    
    %------------------ explicit nonlinear --------------------
    psitilde = rk2gsh1(psimat,zetamat,p,@nlpgsh1); 
    psitilde = latToVec(psitilde);

    omztilde = rk2gsh2(omzmat,psimat,p,@nlpgsh2);
    omztilde = latToVec(omztilde);
    
    %------------------ implicit CN linear ----------------------
    psivec = matdivpsi * psitilde;
    omzvec = matdivomz * omztilde;

    %--------- converting from vector to matrix ----------
    psimat = vecToLat(psivec,Nx,Ny);
    omzmat = vecToLat(omzvec,Nx,Ny);
    
    
    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
    
end

fprintf('Transients done...\n\n');
end

