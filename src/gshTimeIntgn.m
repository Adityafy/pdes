function [psimat, omzmat, zetamat, ulat, vlat] = gshTimeIntgn(p, gshtrfunc)
% [psimat, omzmat, zetamat] = gshTimeIntgn(p, gshtrfunc)
% The time integration function that yields psi, omega_z, zeta and 
% u and v velocities. gshtrfunc (gsh transient time evolution func) 
% yields the initial conditions.
addpath('../src/');

dt = p.dt; % time step size

totimeu = p.totimeu;
tmax = totimeu/dt; % total time steps

Nx = p.Nx;
Ny = p.Ny;
N = Nx*Ny;

seed = p.seed;
rng(seed,"twister");


%% Initial conditions and preallocation
% Running the transients
[psimatic, omzmatic, zetamatic, ~, ~] = gshtrfunc(p);

psimat = zeros(Nx,Ny,1); % scalar field
psivec = zeros(N,1);

psimat(:,:,1) = psimatic;
psivec(:,1) = latToVec(psimat(:,:,1));

zetamat = zeros(Nx,Ny,tmax); % scalar field
ulat = zeros(Nx,Ny,tmax);
vlat = zeros(Nx,Ny,tmax);

zetamat(:,:,1) = zetamatic;

omzmat = zeros(Nx,Ny,tmax);
omzmat(:,:,1) = omzmatic;
omzvec = zeros(N,1);


%% Time integration
run_name = p.run_name;
fprintf(join(['Running ',run_name]));

fprintf('\nCalculating coefficient matrices for implicit calcs...\n');
[matdivpsi, matdivomz] = impMatGSH(p);
toc

fprintf('Running time integration...\n');
for t = 1:tmax    
    %------------------ zeta, iterative ------------------------
    zetamat(:,:,t) = iterativeZetaOmz(omzmat(:,:,t),zetamat(:,:,t),p);
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
    
    
    zetamat(:,:,t+1) = zetamat(:,:,t);
    if rem(t,tmax/10) == 0
        fprintf('This is time step %i, ',t);
        toc
    end
    
end


end

