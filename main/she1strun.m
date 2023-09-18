%% SHE

clear all;
close all;
seed = 1;
rng(seed);

tmax = 5e6;

dt = 0.5e-10;
dspace = 0.5;
dx = dspace;
dy = dspace;

Nx = 20;
Ny = 20;
u = zeros(Nx,Ny,tmax);
u(:,:,1) = rand(Nx,Ny);
% u(1,2,1) = 0.245;

r = 0.9;
%%
% Euler solver with PBC
% I = 1:Nx;
% J = 1:Ny;
% Im1 = circshift(I,1);
% Im2 = circshift(I,2);
% Ip1 = circshift(I,-1);
% Ip2 = circshift(I,-2);
% Jm1 = circshift(J,1);
% Jm2 = circshift(J,2);
% Jp1 = circshift(J,-1);
% Jp2 = circshift(J,-2);
% for t = 1:tmax
%     for i = 1:Nx
%         for j = 1:Ny
%             u(i,j,t+1) = u(i,j,t) + (-2*dt) * ...
%             ( (u(Im1(i),j,t) -2*u(i,j,t) + u(Ip1(i),j,t))/(dx^2) ...
%               + (u(i,Jm1(j),t) -2*u(i,j,t) + u(i,Jp1(j),t))/(dy^2) ) ...
%              - dt * (u(Im2(i),j,t) - 4*u(Im1(i),j,t) + 6*u(i,j,t) - 4*u(Ip1(i),j,t) + u(Ip2(i),j,t))/(dx^4) ...
%              - dt * (u(i,Jm2(j),t) - 4*u(i,Jm1(j),t) + 6*u(i,j,t) - 4*u(i,Jp1(j),t) + u(i,Jp2(j),t))/(dy^4) ...
%              + (r -1)*u(i,j,t) - u(i,j,t)^3;
%         end
%     end
% end
%%
tic
% ufinal = rk4dyn(r,Nx,Ny,dx,dy,u,@she,dt,tmax);
ufinal = foeulerdyn(r,Nx,Ny,dx,dy,u,@she,dt,tmax);
toc
%%
figure;
imagesc(ufinal(:,:,500));
figure;
imagesc(ufinal(:,:,1000));
figure;
imagesc(ufinal(:,:,1500));
figure;
imagesc(ufinal(:,:,tmax));


uVideoFilename = 'uvideo6';
u_video = VideoWriter(uVideoFilename, 'MPEG-4');
% u_video.FrameRate = 60;
open(u_video);
figure();
for t = 1:1000:tmax

    imagesc(ufinal(:,:,t));
    xlabel('i');
    ylabel('j');
    colorbar;
    % clim([-0.6 0.6]);
    % colormap('jet');
    frame = getframe(gcf);
    writeVideo(u_video,frame);
end
close(u_video);


%%
function unp1 = she(r,Nx,Ny,dx,dy,un)
    I = 1:Nx;
    J = 1:Ny;
    Im1 = circshift(I,1);
    Im2 = circshift(I,2);
    Ip1 = circshift(I,-1);
    Ip2 = circshift(I,-2);
    Jm1 = circshift(J,1);
    Jm2 = circshift(J,2);
    Jp1 = circshift(J,-1);
    Jp2 = circshift(J,-2);
    unp1 = zeros(Nx,Ny);
    % putting dt = 1
    for i = 1:Nx
        for j = 1:Ny
            unp1(i,j) = - (un(Im2(i),J(j)) - 4*un(Im1(i),J(j)) + 6*un(I(i),J(j)) ...
                    - 4*un(Ip1(i),J(j)) + un(Ip2(i),J(j)))/(dx^4) ...
                - (un(I(i),Jm2(j)) - 4*un(I(i),Jm1(j)) + 6*un(I(i),J(j)) ...
                    - 4*un(I(i),Jp1(j)) + un(I(i),Jp2(j)))/(dy^4) ...
                    ...
                - 2 * ( un(Ip1(i),Jp1(j)) + un(Ip1(i),Jm1(j)) ...
                    + un(Im1(i),Jm1(j)) + un(Im1(i),Jp1(j)) ) / (dx^2 * dy^2) ...
                - 2 * (-2) * ( un(Ip1(i),J(j)) + un(Im1(i),J(j)) ...
                    + un(I(i),Jm1(j)) + un(I(i),Jp1(j)) ) / (dx^2 * dy^2) ...
                - 2 * (-4) * un(I(i),J(j)) / (dx^2 * dy^2) ...
                ...
                - 2 * ( (un(Im1(i),J(j)) -2*un(I(i),J(j)) + un(Ip1(i),J(j)))/(dx^2) ...
                    + (un(I(i),Jm1(j)) -2*un(I(i),J(j)) + un(I(i),Jp1(j)))/(dy^2) ) ...
                ...
                + (r-1)*un(I(i),J(j)) - un(I(i),J(j))^3;
        end
    end
end

function dynamics = rk4dyn(r,Nx,Ny,dx,dy,dynamics,dynFunc,h,tmax)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
for t = 1:tmax
    % k values
    k1 = dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t));
    k2 = dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t) + (0.5*h)*k1);
    k3 = dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t) + (0.5*h)*k2);
    k4 = dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t) + h*k3);
    % dynamics
    dynamics(:,:,t+1) = dynamics(:,:,t) + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
end
end

function dynamics = foeulerdyn(r,Nx,Ny,dx,dy,dynamics,dynFunc,h,tmax)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
for t = 1:tmax
    % k values
    % k1 = dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t));
    % k2 = dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t) + (0.5*h)*k1);
    % k3 = dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t) + (0.5*h)*k2);
    % k4 = dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t) + h*k3);
    % dynamics
    dynamics(:,:,t+1) = dynamics(:,:,t) + h*dynFunc(r,Nx,Ny,dx,dy, dynamics(:,:,t));
end
end