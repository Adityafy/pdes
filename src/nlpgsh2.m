function np2 = nlpgsh2(~, psi, p)
%NLPGSH2 Computes the nonlinear RHS for the second GSH equation 1 for
%time-step.
%
%   np2 = nlpgsh2(~, psi, p)
%
%   This function calculates the nonlinear right-hand side of the second
%   Generalized Swift-Hohenberg (GSH) equation, used in explicit time
%   integration schemes with periodic boundary conditions.
%
%   Inputs:
%     ~    : Placeholder for unused input (for consistency)
%     psi  : Scalar field psi at current time
%     p    : Struct containing simulation and grid parameters
%
%   Output:
%     np2  : Nonlinear RHS contribution to the evolution of omega_z
%
%   Author: Aditya

    addpath('../src/');  % Ensure required source functions are available

    gm = p.con.gm;        % Nonlinear coupling parameter
    Nx = p.rmesh.Nx;      % Number of grid points in x-direction
    Ny = p.rmesh.Ny;      % Number of grid points in y-direction

    % gmPos = p.gmPos;    % Optional flag for alternate gm placement

    [I, J, ~] = unpackIndexStructPBC(p);  % Periodic indices

    % Initialize output array
    np2 = zeros(Nx, Ny);

    % Loop over spatial domain
    for i = 1:length(I)
        for j = 1:length(J)
            dpsidx        = sfdd(psi, 1, 0, i, j, p);  % ∂ψ/∂x
            dpsidy        = sfdd(psi, 0, 1, i, j, p);  % ∂ψ/∂y
            d3psidx3      = sfdd(psi, 3, 0, i, j, p);  % ∂³ψ/∂x³
            d3psidy3      = sfdd(psi, 0, 3, i, j, p);  % ∂³ψ/∂y³
            d3psidxdy2    = sfdd(psi, 1, 2, i, j, p);  % ∂³ψ/∂x∂y²
            d3psidydx2    = sfdd(psi, 2, 1, i, j, p);  % ∂³ψ/∂y∂x²

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

    %----------------------------
    % Convergence / Stability Checks
    %----------------------------
    if any(isnan(np2), 'all')
        warning('nlpgsh2:NaNDetected', 'NaNs detected in nonlinear RHS np2.');
    end
    if any(isinf(np2), 'all')
        warning('nlpgsh2:InfDetected', 'Infs detected in nonlinear RHS np2.');
    end
    maxVal = max(abs(np2(:)));
    threshold = 1e6;  % Tunable threshold based on the system
    if maxVal > threshold
        warning('nlpgsh2:BlowupDetected', ...
            'Large values detected in nonlinear RHS np2: max = %.3e', maxVal);
    end
end
