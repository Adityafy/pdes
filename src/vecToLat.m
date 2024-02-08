function xlat = vecToLat(x,Nx,Ny)
% xlat = vecToLat(x,Nx,Ny)
% Converts an Nx*Ny length vector into a (Nx X Ny) size matrix.
    if length(x) ~= Nx*Ny
        msg = "Output lattice dimensions do not comply " + ...
            "with the length of input vector";
        error(msg);
    end
    xlat = zeros(Nx,Ny);
    for i = 1:Nx*Ny
        q = fix(i/Ny);
        r = rem(i,Ny);
        k1 = q+1;
        k2 = r;
        if r == 0
            k1 = q;
            k2 = Ny;
        end
        % xlat(k1,k2) = x(i,i);
        xlat(k1,k2) = x(i);
    end
end