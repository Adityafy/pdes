function np1 = nlpgsh1(psit,zetat,p)
% np1 = nlpgsh1(psit,zetat,p)
% Calculates the RHS of the nonlinear part of the first GSH equation for 
% explicit time marching methods.
% (periodic boundary conditions)

    addpath('../src/');
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

    % c = 1/(4*dx*dy);
    np1 = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)

            %----previous----
            % dzetadx = (zetat(Ip1(i),J(j)) - zetat(Im1(i),J(j)))/(2*dx);
            % dzetady = (zetat(I(i),Jp1(j)) - zetat(I(i),Jm1(j)))/(2*dy);
            % dpsidx = (psit(Ip1(i),J(j)) - psit(Im1(i),J(j)))/(2*dx);
            % dpsidy = (psit(I(i),Jp1(j)) - psit(I(i),Jm1(j)))/(2*dy);
            % psi = psit(I(i),J(j));

            %----new----
            dzetadx = sfdd(zetat,1,0,i,j,p);
            dzetady = sfdd(zetat,0,1,i,j,p);
            dpsidx = sfdd(psit,1,0,i,j,p);
            dpsidy = sfdd(psit,0,1,i,j,p);
            psi = psit(I(i),J(j));

            % if gmPos == 1
            %     np1(i,j) = - 3*psi^3 - gm* ( dzetady * dpsidx - dzetadx * dpsidy );
            % elseif gmPos == 2
                np1(i,j) = - psi^3 - dzetady * dpsidx + dzetadx * dpsidy ;
            % else
            %     error('gmPos must be either 1 or 2');
            % end
        end
    end
end