function idx = makeIndexStructPBC(Nx, Ny)
%MAKEINDEXSTRUCTPBC Creates index arrays for 2D grids with periodic boundary conditions.
%
%   idx = MAKEINDEXSTRUCTPBC(Nx, Ny) returns a struct with fields:
%   I, J, Ip1, Im1, Ip2, Im2, Jp1, Jm1, Jp2, Jm2
%   representing current and shifted indices with periodic wrapping.

    idx.I = 1:Nx;
    idx.J = 1:Ny;

    idx.Ip1 = circshift(idx.I, -1);
    idx.Im1 = circshift(idx.I, 1);
    idx.Ip2 = circshift(idx.I, -2);
    idx.Im2 = circshift(idx.I, 2);

    idx.Jp1 = circshift(idx.J, -1);
    idx.Jm1 = circshift(idx.J, 1);
    idx.Jp2 = circshift(idx.J, -2);
    idx.Jm2 = circshift(idx.J, 2);
end
