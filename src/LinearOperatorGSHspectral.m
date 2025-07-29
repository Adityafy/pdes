function [L1,L2] = LinearOperatorGSHspectral(p,Kx,Ky)
% [L1,L2] = LinearOperatorGSHspectral(p,Kx,Ky)
% L1 and L2 are the linear operator matrices for GSH equation 1 and 2
% respectively.
epsilon = p.con.epsilon;
sig = p.con.sig;
c = p.con.c;
L1 = epsilon - 1 - Kx.^4 - 2*(Kx.^2).*(Ky.^2) - Ky.^4 + 2*Kx.^2 + 2*Ky.^2 ;
L2 = - sig * (Kx.^2 + Ky.^2) - sig * c^2;
end