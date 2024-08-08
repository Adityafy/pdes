function [matdivpsi, matdivomz] = impMatGSH(p,h)
% [matdivpsi, matdivomz] = impMatGSH(p)
% Calculates coefficient matrices for implicit calculations in GSH
% equations 1 and 2.
    


    %---------------from parameters-------------
    eps = p.eps;
    sig = p.sig;
    c = p.c;
    
    Nx = p.Nx;
    Ny = p.Ny;
    N = Nx*Ny;
    % h = p.dt;
    delta = p.dx;
    dx = p.dx;
    dy = p.dy;

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

    %----------coefficients--------------
    % zeta
    % z1 = 1/dx^2;
    % z2 = -2 * ((1/dx^2) + (1/dy^2));
    % z3 = 1/dy^2;
    
    
    a0 = - (eps-1)*h/2 + 10*h/delta^4 - 4*h/delta^2;
    
    a1 = 1 + a0;
    a2 = -4*h/delta^4 + h/delta^2;
    a3 = h/delta^4;
    a4 = h/(2*delta^4);
    
    b1 = 1 - a0;
    b2 = -a2;
    b3 = -a3;
    b4 = -a4;
    
    % omz
    m0 = sig*h/(dx^2) + sig*h/(dy^2) + h*sig*(c^2)/2;
    m1 = 1 + m0;
    m2 = - sig*h/(2*dx^2);
    m3 = - sig*h/(2*dy^2);
    m4 = 1 - m0;
    m5 = - m2;
    m6 = - m3;
    
    % TS 2
    e0 = -m0;
    e1 = 1+e0;
    e2 = -m2;
    e3 = -m3;
    
    
    
    %---------------preallocation-----------------
    % z = zeros(Nx,Ny);
    % Z = zeros(N,N);
    
    a = zeros(Nx,Ny);
    b = zeros(Nx,Ny);
    A = zeros(N,N);
    B = zeros(N,N);
    
    ml = zeros(Nx,Ny);
    mr = zeros(Nx,Ny);
    Ml = zeros(N,N);
    Mr = zeros(N,N);
    

    %------------calculation-----------------
    n = 1;
    for i = 1:length(I)
        for j = 1:length(J)
            % z(I(i),J(j)) = z2;
            % z(I(i),Jp1(j)) = z3;
            % z(I(i),Jm1(j)) = z3;
            % z(Ip1(i),J(j)) = z1;
            % z(Im1(i),J(j)) = z1;
            % Z(n,:) = latToVec(z)';
    
            a(I(i),J(j)) = a1;
    
            a(I(i),Jp1(j)) = a2;
            a(I(i),Jm1(j)) = a2;
            a(Ip1(i),J(j)) = a2;
            a(Im1(i),J(j)) = a2;
    
            a(Ip1(i),Jp1(j)) = a3;
            a(Ip1(i),Jm1(j)) = a3;
            a(Im1(i),Jp1(j)) = a3;
            a(Im1(i),Jm1(j)) = a3;
    
            a(Im2(i),J(j)) = a4;
            a(I(i),Jp2(j)) = a4;
            a(I(i),Jm2(j)) = a4;
            a(Ip2(i),J(j)) = a4;
    
            A(n,:) = latToVec(a)';
    
            b = -1*a;
    
            b(I(i),J(j)) = b1;
    
            B(n,:) = latToVec(b)';
            
            ml(I(i),J(j)) = m1;
            ml(Im1(i),J(j)) = m2;
            ml(Ip1(i),J(j)) = m2;
            ml(I(i),Jm1(j)) = m3;
            ml(I(i),Jp1(j)) = m3;
    
            Ml(n,:) = latToVec(ml)';
    
            mr = - ml;
    
            mr(I(i),J(j)) = 1 - m0;
    
            Mr(n,:) = latToVec(mr)';
    
            n = n+1;
            % z = zeros(Nx,Ny);
            a = zeros(Nx,Ny);
            b = zeros(Nx,Ny);
            ml = zeros(Nx,Ny);
            mr = zeros(Nx,Ny);
        end
    end
    
    % Z = sparse(Z);
    
    % Lz = tril(Z);
    % Uz = triu(Z,1);
    % Lzinv = inv(Lz);
    
    matdivpsi = A\B;
    matdivomz = Ml\Mr;

end