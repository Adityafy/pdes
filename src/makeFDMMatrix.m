function fdA = makeFDMMatrix(p)
% makeFDMMatrix Construct finite-difference Laplacian matrix with PBC
%   fdA = makeFDMMatrix(p) returns the sparse matrix representing the 2D
%   Laplacian with periodic boundaries.
%
%   Input:
%       p: struct with fields
%           - Nx, Ny : number of grid points in x and y directions
%           - dx     : grid spacing (assumed same in both directions)
%           - N      : total number of grid points (Nx * Ny)
%
%   Output:
%       fdA: sparse matrix of size (N x N)

    % Preallocate matrix
    fdA = zeros(p.rmesh.N);
    delsq = p.rmesh.dx^2;

    % Index helpers
    I = 1:p.rmesh.Nx;
    J = 1:p.rmesh.Ny;
    Ip1 = circshift(I, -1);
    Im1 = circshift(I, 1);
    Jp1 = circshift(J, -1);
    Jm1 = circshift(J, 1);

    % Matrix construction
    spotter = zeros(p.Nx, p.Ny);
    iters = 1;
    for i = 1:p.Nx
        for j = 1:p.Ny
            spotter(i, j) =  4 / delsq;
            spotter(Ip1(i), j) = -1 / delsq;
            spotter(Im1(i), j) = -1 / delsq;
            spotter(i, Jp1(j)) = -1 / delsq;
            spotter(i, Jm1(j)) = -1 / delsq;

            fdA(iters, :) = reshape(spotter', [], 1)';  % 2D to row vector
            iters = iters + 1;

            spotter(:) = 0;  % Reset for next iteration
        end
    end

    % Return as sparse matrix
    fdA = sparse(fdA);
end
