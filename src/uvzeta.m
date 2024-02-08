function [u,v] = uvzeta(zeta,p)
% [u,v] = uvzeta(zeta,p)
% Part of the gsh calculations.
% Calculates u and v velocities from zeta.
% (periodic boundary conditions)
    dx = p.dx;
    dy = p.dy;
    Nx = p.Nx;
    Ny = p.Ny;

    I = p.I;
    J = p.J;
    Ip1 = p.Ip1;
    Im1 = p.Im1;
    Jp1 = p.Jp1;
    Jm1 = p.Jm1;

    u = zeros(Nx,Ny);
    v = zeros(Nx,Ny);
    for i = 1:length(I)
        for j = 1:length(J)
            dzetadx = (zeta(Ip1(i),J(j)) - zeta(Im1(i),J(j)))/(2*dx);
            dzetady = (zeta(I(i),Jp1(j)) - zeta(I(i),Jm1(j)))/(2*dy);
            u(i,j) = dzetady;
            v(i,j) = -dzetadx;
        end
    end
end