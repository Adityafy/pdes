function zeta = iterativeZetaOmz(omzmat, zetamat, p)
% zeta = iterativeZetaOmz(omzmat, zetamat, p)
% An iterative poisson solver for
% Omega_z = - Laplacian(Zeta)
% with periodic boundary conditions.
% mind the negative sign if used for other kinds of Poisson solvers
% Number of iterations is fixed for now (= 800)
    delta = p.dx;

    I = p.I;
    J = p.J;
    Ip1 = p.Ip1;
    Im1 = p.Im1;
    Jp1 = p.Jp1;
    Jm1 = p.Jm1;

    zeta = zetamat;
    iterations = 800;
    difference = zeros(1,iterations);

    for k = 1:iterations
        for i = 1:length(I)
            for j = 1:length(J)
                zeta(i,j) = ((delta^2)/4) * omzmat(I(i),J(j)) ...
                    + (1/4) * ( zetamat(Im1(i),J(j)) + zetamat(Ip1(i),J(j)) ...
                    + zetamat(I(i),Jm1(j)) + zetamat(I(i),Jp1(j)) );
            end
        end

        lapl = zeros(size(zeta));
        for i = 1:length(I)
            for j = 1:length(J)
                lapl(i,j) = (1/(delta^2)) * ...
                    ( zeta(Im1(i),J(j)) + zeta(Ip1(i),J(j)) ...
                    - 4 * zeta(I(i),J(j)) + ...
                    + zeta(I(i),Jm1(j)) + zeta(I(i),Jp1(j)) );
            end
        end
        b = lapl - (-omzmat);
        difference(k) = max(max(abs(b)));
        zetamat = zeta;
    end
end