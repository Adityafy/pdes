function [dpsi1,domz1,fpv] = gshfpvic(p)

addpath('../src/');
Nx = p.Nx;
Ny = p.Ny;
N = p.N;

dpsi1 = zeros(Nx,Ny);

domz1 = zeros(Nx,Ny); 

if p.fpvictype < 0 % this is default
    % dpsi1(1,1,1) = 0.001;
    dpsi1 = -0.05 + (0.05-(-0.05)) * rand(Nx,Ny);
    rng(p.seed+2,"twister"); % setting different IC for omz
    domz1 = -0.05 + (0.05-(-0.05)) * rand(Nx,Ny);
    rng(p.seed, "twister"); % putting back default seed
    fpv = [latToVec(dpsi1(:,:,1)); latToVec(domz1(:,:,1))];
    fpv(:,1) = fpv(:,1)./ norm(fpv(:,1));
elseif p.fpvictype == 0 % all zeros
    dpsi1(:,:,1) = zeros(Nx,Ny);
    domz1(:,:,1) = zeros(Nx,Ny);
    fpv = [latToVec(dpsi1(:,:,1)); latToVec(domz1(:,:,1))];
elseif p.fpvictype > 0 % all same value
    dpsi1(:,:,1) = p.fpvictype * ones(Nx,Ny);
    domz1(:,:,1) = p.fpvictype * ones(Nx,Ny);
    fpv = [latToVec(dpsi1(:,:,1)); latToVec(domz1(:,:,1))];
else
    error('fpv ic type should be a scalar.')
end

dpsi1(:,:,1) = vecToLat(fpv(1:N,1),Nx,Ny);
domz1(:,:,1) = vecToLat(fpv(N+1:2*N,1),Nx,Ny);

end