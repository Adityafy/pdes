function np1 = nlpgsh1(psit, zetat, p)
%NLPGSH1 Computes the nonlinear RHS for the first GSH equation for one time
%step.
%
%   np1 = nlpgsh1(psit, zetat, p)
%
%   This function calculates the nonlinear right-hand side of the first
%   Generalized Swift-Hohenberg (GSH) equation, used in explicit time
%   marching schemes with periodic boundary conditions.
%
%   Inputs:
%     psit   : Scalar field psi at current time
%     zetat  : Scalar field zeta at current time
%     p      : Struct containing grid and simulation parameters
%
%   Output:
%     np1    : Nonlinear RHS contribution to the evolution of psi
%
%   Author: Aditya

    addpath('../src/');  % Ensure required source functions (e.g. sfdd) are available

    Nx = p.rmesh.Nx;     % Number of grid points in x-direction
    Ny = p.rmesh.Ny;     % Number of grid points in y-direction

    % gmPos = p.con.gmPos;   % (optional: kept for reference)

    [I, J, ~] = unpackIndexStructPBC(p);  % Periodic grid indices

    % Initialize output array
    np1 = zeros(Nx, Ny);

    % Loop over all spatial points
    for i = 1:length(I)
        for j = 1:length(J)
            %----new----
            dzetadx = sfdd(zetat, 1, 0, i, j, p);  % ∂zeta/∂x
            dzetady = sfdd(zetat, 0, 1, i, j, p);  % ∂zeta/∂y
            dpsidx  = sfdd(psit,  1, 0, i, j, p);  % ∂psi/∂x
            dpsidy  = sfdd(psit,  0, 1, i, j, p);  % ∂psi/∂y
            psi     = psit(I(i), J(j));            % Current psi value

            % if gmPos == 1
            %     np1(i,j) = - 3*psi^3 - gm* ( dzetady * dpsidx - dzetadx * dpsidy );
            % elseif gmPos == 2
            np1(i,j) = - psi^3 - dzetady * dpsidx + dzetadx * dpsidy ;
            % else
            %     error('gmPos must be either 1 or 2');
            % end
        end
    end

    %----------------------------
    % Convergence / Stability Checks
    %----------------------------
    if any(isnan(np1), 'all')
        warning('nlpgsh1:NaNDetected', 'NaNs detected in nonlinear RHS np1.');
    end
    if any(isinf(np1), 'all')
        warning('nlpgsh1:InfDetected', 'Infs detected in nonlinear RHS np1.');
    end
    maxVal = max(abs(np1(:)));
    threshold = 1e6;  % Adjust as needed based on expected dynamic range
    if maxVal > threshold
        warning('nlpgsh1:BlowupDetected', ...
            'Large values detected in nonlinear RHS np1: max = %.3e', maxVal);
    end
end
