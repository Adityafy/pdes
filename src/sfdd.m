function derivtv = sfdd(v,xd,yd,i,j,p)
% derivtv = sfdd(v,xd,yd,i,j,p)
% Yields the spatial finite difference derivative (sfdd) for GSH calcs.
% 2nd order accurate central differences.
% Periodic Boundary Conditions!!!
% Not all derivatives are possible.
% Example: d3v/dxdy2(i,j) = sfdd(v,1,2,i,j,p)
dx = p.rmesh.dx;
dy = p.rmesh.dy;

[I, J, Ip1, Im1, Jp1, Jm1, Ip2, Im2, Jp2, Jm2] = unpackIndexStructPBC(p);

% for i = 1:length(I)
%     for j = 1:length(J)

if xd == 1 && yd == 0 % dv/dx
    derivtv = v(Ip1(i),J(j))/(2*dx) - v(Im1(i),J(j))/(2*dx);
elseif xd == 0 && yd == 1 % dv/dy
    derivtv = v(I(i),Jp1(j))/(2*dy) - v(I(i),Jm1(j))/(2*dy);
elseif xd == 2 && yd == 0
    derivtv = (1/(dx^2)) * (v(Im1(i),J(j)) - 2 * v(I(i),J(j)) ...
        + v(Ip1(i),J(j)));
elseif xd == 0 && yd == 2
    derivtv = (1/(dy^2)) * (v(I(i),Jm1(j)) - 2 * v(I(i),J(j)) ...
        + v(I(i),Jp1(j)));
elseif xd == 3 && yd == 0 % d3v/dx3
    derivtv = (1/(2*dx^3)) * (-v(Im2(i),J(j)) + 2 * v(Im1(i),J(j)) ...
        - 2 * v(Ip1(i),J(j)) + v(Ip2(i),J(j)));
elseif xd == 0 && yd == 3 % d3v/dy3
    derivtv = (1/(2*dy^3)) * (-v(I(i),Jm2(j)) + 2 * v(I(i),Jm1(j)) ...
        - 2 * v(I(i),Jp1(j)) + v(I(i),Jp2(j)));
elseif xd == 1 && yd == 2 % d3v/dxdy2
    derivtv = (1/(2*dx*dy^2)) * (v(Ip1(i),Jm1(j)) - v(Im1(i),Jm1(j)) ...
        - 2 * v(Ip1(i),J(j)) + 2*v(Im1(i),J(j)) ...
        + v(Ip1(i),Jp1(j)) - v(Im1(i),Jp1(j)));
elseif xd == 2 && yd == 1 % d3v/dydx2
    derivtv = (1/(2*dy*dx^2)) * (v(Im1(i),Jp1(j)) - v(Im1(i),Jm1(j)) ...
        - 2 * v(I(i),Jp1(j)) + 2 * v(I(i),Jm1(j)) ...
        + v(Ip1(i),Jp1(j)) - v(Ip1(i),Jm1(j)));
elseif xd == 4 && yd == 0 %d4v/dx4
    derivtv = (1/(dx^4)) * (v(Im2(i),J(j)) - 4 * v(Im1(i),J(j)) ...
        + 6*v(I(i),J(j)) - 4 * v(Ip1(i),J(j)) + v(Ip2(i),J(j)));
elseif xd == 0 && yd == 4 % d4v/dy4
    derivtv = (1/(dy^4)) * (v(I(i),Jm2(j)) - 4 * v(I(i),Jm1(j)) ...
        + 6*v(I(i),J(j)) - 4 * v(I(i),Jp1(j)) + v(I(i),Jp2(j)));
elseif xd == 2 && yd == 2 % d4v/dx2dy2
    derivtv = (1/(dx^2*dy^2)) * ( v(Ip1(i),Jp1(j)) + v(Ip1(i),Jm1(j))  ...
        + v(Im1(i),Jp1(j)) + v(Im1(i),Jm1(j))) ...
        - (2/(dx^2*dy^2)) * ( v(Ip1(i),J(j)) + v(Im1(i),J(j)) ...
        + v(I(i),Jm1(j)) +  v(I(i),Jp1(j)) ) ...
        + (4/(dx^2*dy^2)) * ( v(I(i),J(j)) );
else
    error(['Requested spatial finite difference ' ...
        'derivative is not available!']);
end
%     end
% end
end