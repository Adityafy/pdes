function np2 = nlpgsh2(~,psi,p)
% np1 = nlpgsh2(psit,zetat,p)
% Calculates the RHS of the nonlinear part of the second GSH equation for 
% explicit time marching methods.
% (periodic boundary conditions)
    
    addpath('../src/');
    gm = p.gm;
    dx = p.dx;
    dy = p.dy;
    Nx = p.Nx;
    Ny = p.Ny;
    gmPos = p.gmPos;

    I = p.I;
    J = p.J;
    Ip1 = p.Ip1;
    Im1 = p.Im1;
    Jp1 = p.Jp1;
    Jm1 = p.Jm1;
    
    Ip2 = p.Ip2;
    Im2 = p.Im2;
    Jp2 = p.Jp2;
    Jm2 = p.Jm2;

    np2 = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)

            %----previous----
            % dpsidx = psi(Ip1(i),J(j))/(2*dx) - psi(Im1(i),J(j))/(2*dx);
            % dpsidy = psi(I(i),Jp1(j))/(2*dy) - psi(I(i),Jm1(j))/(2*dy);
            % d3psidx3 = (1/(2*dx^3)) * (-psi(Im2(i),J(j)) + 2 * psi(Im1(i),J(j)) ...
            %         - 2 * psi(Ip1(i),J(j)) + psi(Ip2(i),J(j)));
            % d3psidy3 = (1/(2*dy^3)) * (-psi(I(i),Jm2(j)) + 2 * psi(I(i),Jm1(j)) ...
            %         - 2 * psi(I(i),Jp1(j)) + psi(I(i),Jp2(j)));
            % d3psidxdy2 = (1/(2*dx*dy^2)) * (psi(Ip1(i),Jm1(j)) - psi(Im1(i),Jm1(j)) ...
            %         - 2 * psi(Ip1(i),J(j)) + 2*psi(Im1(i),J(j)) ...
            %         + psi(Ip1(i),Jp1(j)) - psi(Im1(i),Jp1(j)));
            % d3psidydx2 = (1/(2*dy*dx^2)) * (psi(Im1(i),Jp1(j)) - psi(Im1(i),Jm1(j)) ...
            %         - 2 * psi(I(i),Jp1(j)) + 2 * psi(I(i),Jm1(j)) ...
            %         + psi(Ip1(i),Jp1(j)) - psi(Ip1(i),Jm1(j)));

            %----new----
            dpsidx = sfdd(psi,1,0,i,j,p);
            dpsidy = sfdd(psi,0,1,i,j,p);
            d3psidx3 = sfdd(psi,3,0,i,j,p);
            d3psidy3 = sfdd(psi,0,3,i,j,p);
            d3psidxdy2 = sfdd(psi,1,2,i,j,p);
            d3psidydx2 = sfdd(psi,2,1,i,j,p);

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
end