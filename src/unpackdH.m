function [dpsimat, dpsihmat, domzmat, domzhmat, ...
                    dzetamat, dH, dHmag] = unpackdH(p,dH)

nv = p.ts.nv;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;
N = p.rmesh.N;
nmax = p.sim.nmax;
dHmag = zeros(nv,nmax); % magnitude of the perturbation vectors
dpsimat = zeros(Nx,Ny,nv);
domzmat = zeros(Nx,Ny,nv);
dzetamat = zeros(Nx,Ny,nv);
dpsihmat = zeros(size(dpsimat));
domzhmat = zeros(size(domzmat));

% same as fpv fd si
rng(p.sim.seed,"twister");
% randrange = 0.01;
for k = 1:nv
    rng(p.sim.seed+k,"twister");
    % dH(:,k,1) = -randrange + (randrange-(-randrange)) * rand(2*N,1,1);
    dH(:,k,1) = dH(:,k,1)./norm(dH(:,k,1));
    dpsimat(:,:,k) = reshape(dH(1:N,k,1),Nx,Ny)';
    domzmat(:,:,k) = reshape(dH(N+1:2*N,k,1),Nx,Ny)';
    dzetamat(:,:,k) = zetaGSHspectral(p,domzmat(:,:,k));
    dHmag(k,1) = norm(dH(:,k,1));
end

rng(p.sim.seed,"twister"); % set rng back to seed

for k = 1:nv
    dpsihmat(:,:,k) = fft2(dpsimat(:,:,k));
    domzhmat(:,:,k) = fft2(domzmat(:,:,k));
end

laminst = zeros(nv,1); % lambda(k) instantaneous


end