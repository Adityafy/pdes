% script for the equivalence of reshape (matlab) and lattovec (created by Aditya)

addpath('../src/');

Nx = 4;
Ny = 4;
A = magic(Nx);
Avec = latToVec(A);

fprintf("Is latToVec(A) equal to reshape(A',[],1)?")
isequal(latToVec(A),reshape(A',[],1))
fprintf("\nIs vecToLat(Avec,Nx,Ny) equal to reshape(Avec,Nx,Ny)'?")
isequal(vecToLat(Avec,Nx,Ny),reshape(Avec,Nx,Ny)')