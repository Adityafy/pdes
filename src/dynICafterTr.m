function [psimat,psi,omzmat,omz,zetamat,zeta] = dynICafterTr(dynAddress)
% [psimat,psi,omzmat,omz,zetamat,zeta] = dynICafterTr(dynAddress)

fprintf('\n Loading initial conditions.\n');
load(dynAddress,'psiICpostTr','omzICpostTr', 'zetaICpostTr', ...
    'uICpostTr','vICpostTr');
toc;
psimat = psiICpostTr; % assigning end of transients to the 
                             % beginning of dynamics after transients
zetamat = zetaICpostTr;
omzmat = omzICpostTr;

%%%======= Initializing variables for time marching ========
psi(:,:,1) = psimat; % 3D mat that's going to be appended and saved
% psihmat = fft2(psimat); % assigning psihat for fft of psi

zeta(:,:,1) = zetamat;
umat = uICpostTr;
vmat = vICpostTr;
% u(:,:,1) = umat;
% v(:,:,1) = vmat;
% zetahmat = fft2(zetamat);

% omztr = -fdlaplacian(zetatr, Nx,Ny,dx); % if fd laplacian is to be used
% omzmat = (Kx.^2+Ky.^2).*(zetatrmat);
omz(:,:,1) = omzmat;
% omzhmat = fft2(omzmat);
end