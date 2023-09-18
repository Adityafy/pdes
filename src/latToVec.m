%% Function to convert lattice matrix into diagonal matrix
function xvec = latToVec(x)
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
