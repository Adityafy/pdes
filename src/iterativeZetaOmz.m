function zeta = iterativeZetaOmz(omzmat, zetamat, p)
%ITERATIVEZETAOMZ Iterative Poisson solver for vorticity-streamfunction relation.
%
%   zeta = iterativeZetaOmz(omzmat, zetamat, p)
%
%   Solves the Poisson equation:
%       Omega_z = -Laplacian(Zeta)
%   with periodic boundary conditions using a fixed number of iterations.
%
%   Inputs:
%     omzmat   : Vorticity field (Omega_z)
%     zetamat  : Initial guess for streamfunction (Zeta)
%     p        : Struct containing grid spacing and index arrays
%
%   Output:
%     zeta     : Updated streamfunction field after iterations
%
%   Notes:
%   - Uses finite difference Laplacian on a periodic grid.
%   - Negative sign in the Poisson equation is important for conventions.
%   - Number of iterations is currently fixed to 800.
%
%   Author: Aditya

    delta = p.rmesh.dx;  % Grid spacing

    % Unpack periodic index arrays
    [I, J, Ip1, Im1, Jp1, Jm1, ~, ~, ~, ~] = unpackIndexStructPBC(p);

    zeta = zetamat;              % Initialize zeta
    iterations = 800;            % Fixed number of iterations
    difference = zeros(1, iterations);  % Optional convergence monitor (not used further)

    for k = 1:iterations
        % Gauss-Seidel update of zeta using 5-point stencil
        for i = 1:length(I)
            for j = 1:length(J)
                zeta(i,j) = (delta^2 / 4) * omzmat(I(i), J(j)) + (1/4) * ( ...
                    zetamat(Im1(i), J(j)) + zetamat(Ip1(i), J(j)) + ...
                    zetamat(I(i), Jm1(j)) + zetamat(I(i), Jp1(j)) );
            end
        end

        % Compute Laplacian of updated zeta
        lapl = zeros(size(zeta));
        for i = 1:length(I)
            for j = 1:length(J)
                lapl(i,j) = (1 / delta^2) * ( ...
                    zeta(Im1(i), J(j)) + zeta(Ip1(i), J(j)) - ...
                    4 * zeta(I(i), J(j)) + ...
                    zeta(I(i), Jm1(j)) + zeta(I(i), Jp1(j)) );
            end
        end

        % Residual (not used to terminate, just for monitoring)
        b = lapl + omzmat;
        difference(k) = max(max(abs(b)));

        % Update zetamat for next iteration
        zetamat = zeta;
    end
end
