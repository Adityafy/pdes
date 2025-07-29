function [dpsimat, dpsihmat, domzmat, domzhmat, ...
                    dzetamat, dH, laminst, dHmag] = tsics(p)
nv = p.ts.nv;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;
Kx = p.smesh.Kx;
Ky = p.smesh.Ky;
N = p.rmesh.N;
nmax = p.sim.nmax;

% same as fpv fd si
rng(p.sim.seed,"twister");
randrange = 0.05;
for k = 1:nv
    rng(p.sim.seed+k,"twister");
    dpsimat(:,:,k) = -randrange + (randrange-(-randrange)) * rand(Nx,Ny);
    dpsimat(:,:,k) = dpsimat(:,:,k)./norm(dpsimat(:,:,k));
    rng(p.sim.seed+2+k,"twister"); % setting different IC for omz
    domzmat(:,:,k) = -randrange + (randrange-(-randrange)) * rand(Nx,Ny);
    domzmat(:,:,k) = domzmat(:,:,k)./norm(domzmat(:,:,k));
    dzetamat(:,:,1) = zetaGSHspectral(p,domzmat(:,:,k));
end
rng(p.sim.seed,"twister");

% original process
% rng(seed + 2, "twister");
% randrange = 0.05;
% dpsimat = -randrange + (randrange + randrange) * rand(Nx,Ny,nv);
% dzetamat = -randrange + (randrange + randrange) * rand(Nx,Ny,nv);
% domzmat = (Kx.^2+Ky.^2) .* dzetamat;

% dpsivec = latToVec(dpsimat);
% domzvec = latToVec(domzmat);

dpsihmat = zeros(size(dpsimat));
domzhmat = zeros(size(domzmat));
dH = zeros(2*N,nv,1); % columns of dH are perturbation vectors
% dH = 0.001*orth(dH); % make them orthonormal and small

for k = 1:nv
    dpsihmat(:,:,k) = fft2(dpsimat(:,:,k));
    domzmat(:,:,k) = (Kx.^2+Ky.^2).*(dzetamat(:,:,k));
    domzhmat(:,:,k) = fft2(domzmat(:,:,k));
    % dpsivec(:,k) = latToVec(dpsimat(:,:,k));
    % domzvec(:,k) = latToVec(domzmat(:,:,k));
    dH(:,k,1) = [reshape(dpsimat(:,:,k)',[],1); ...
                 reshape(domzmat(:,:,k)',[],1)];
    % dHhat(:,k,1)= [reshape(dpsihmat(:,:,k)',[],1); ...
    %              reshape(domzhmat(:,:,k)',[],1)];
end
% dH(:,:,1) = [dpsivec; domzvec]; % initializing perturbation vectors

laminst = zeros(nv,1); % lambda(k) instantaneous

dHmag = zeros(nv,nmax); % magnitude of the perturbation vectors
end