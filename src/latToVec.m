function xvec = latToVec(x)
% xvec = latToVec(x)
% Converts a (Nx X Ny) matrix into a column vector of size (Nx*Ny X 1)
    [Nx,Ny] = size(x);
    xvec = zeros(Nx*Ny,1);
    for i = 1:Nx
        for j = 1:Ny
            k = (i-1)*Ny+j;
            xvec(k) = x(i,j);
        end
    end
    % xvec = xvec';
    % xdiag = diag(xvec);
end