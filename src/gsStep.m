function [dpsimat, domzmat, dpsihmat, domzhmat, laminst] = ...
    gsStep(p, n, dpsimat, domzmat, dpsihmat, domzhmat, dH, dHmag)

% Inputs:
%   p         : struct with fields ts.nv, ts.nnorm, grid.Nx, grid.Ny
%   n         : current time step
%   dpsimat   : [Nx x Ny x nv] tangent perturbations in ψ field (real space)
%   domzmat   : [Nx x Ny x nv] tangent perturbations in ω field (real space)
%   dpsihmat  : Fourier transformed ψ perturbations
%   domzhmat  : Fourier transformed ω perturbations
%   dH        : [2*N x nv] tangent vectors stacked (N = Nx*Ny)
%   dHmag     : [nv x ntotal] magnitudes of perturbations
%   laminst   : [nv x numRenorms] instantaneous LEs
%   dH1       : [2*N x numRenorms] first CLV direction stored
%   tN        : renormalization interval (time)

% Outputs:
%   updated versions of all inputs

nv = p.ts.nv;
nnorm = p.ts.nnorm;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;
N = p.rmesh.N;
tN = p.ts.tN;

% for k = 1:nv
%     dH(:,k) = [reshape(real(dpsimat(:,:,k))', [], 1);
%                reshape(real(domzmat(:,:,k))', [], 1)];
%     dHmag(k,n) = norm(dH(:,k));
% end

% --------- Renormalization and LLE Calculation ----------
if rem(n, nnorm) == 0
    [Q, R] = qr(dH);
    dH = Q(:, 1:nv);
    dH1(:, 1) = dH(:,1); % Save first CLV

    for k = 1:nv
        laminst(k, 1) = (1/tN) * log(abs(R(k,k)));
    end

    for k = 1:nv
        dpsimat(:,:,k)  = reshape(dH(1:N,k), Nx, Ny)';
        domzmat(:,:,k)  = reshape(dH(N+1:end,k), Nx, Ny)';
        dpsihmat(:,:,k) = fft2(dpsimat(:,:,k));
        domzhmat(:,:,k) = fft2(domzmat(:,:,k));
        dHmag(k,1) = norm(dH(:,k)); % redundant but retained
    end
end
end
