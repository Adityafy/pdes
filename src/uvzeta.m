function [u, v] = uvzeta(zeta, p)
%UVZETA Computes velocity components u and v from streamfunction zeta.
%
%   [u, v] = uvzeta(zeta, p)
%
%   This function calculates the velocity field (u, v) from the
%   streamfunction zeta using centered finite differences and
%   periodic boundary conditions:
%
%       u =  d(zeta)/dy
%       v = -d(zeta)/dx
%
%   Inputs:
%     zeta : Streamfunction field (2D matrix)
%     p    : Parameter struct containing mesh/grid information
%
%   Outputs:
%     u, v : Velocity components on the same grid as zeta
%
%   Author: Aditya

    % Extract grid spacing and grid size
    dx = p.rmesh.dx;
    dy = p.rmesh.dy;
    Nx = p.rmesh.Nx;
    Ny = p.rmesh.Ny;

    % Unpack periodic indices for finite difference stencils
    [I, J, Ip1, Im1, Jp1, Jm1, ~, ~, ~, ~] = unpackIndexStructPBC(p);

    % Initialize velocity fields
    u = zeros(Nx, Ny);
    v = zeros(Nx, Ny);

    % Loop over interior grid points using periodic indices
    for i = 1:length(I)
        for j = 1:length(J)
            % Centered differences
            dzetadx = (zeta(Ip1(i), J(j)) - zeta(Im1(i), J(j))) / (2 * dx);
            dzetady = (zeta(I(i), Jp1(j)) - zeta(I(i), Jm1(j))) / (2 * dy);

            % Compute velocities from zeta
            u(i, j) =  dzetady;
            v(i, j) = -dzetadx;
        end
    end
end
