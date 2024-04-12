

%% dynamics video

altframes = 2;

    dynvideoname = 'gshDyn';
    dynVideoFilename = join([dynvideoname p.run_name]);
    lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
    %lat_dyn_video.FrameRate = 30;
    open(lat_dyn_video);
    [X,Y] = meshgrid(1:p.Nx,1:p.Ny);
    figure();
    % hold on;
    for t = 1:1:length(psi(1,1,:))-1
        % imagesc(psilat(:,:,t));
        plot1 = contourf(psi(:,:,t),'levels',0.1, 'Linecolor', 'none');
        set(gca,'YDir','normal');
        hold on;
        % imagesc(zetalat(:,:,t));
        colorbar;
        colormap jet
        clim([-0.8 0.8]);
        plot2 = quiver(X,Y,v(:,:,t),u(:,:,t),2,'black');
        xlim([1 p.Nx]);
        ylim([1 p.Ny]);
        frame = getframe(gcf);
        writeVideo(lat_dyn_video,frame);
        clear plot1 plot2;
        hold off;
    end
    hold off;
    close(lat_dyn_video);


%% dpsi1 video with rolls
altframes = 2;
dynvideoname = 'gshrollFPVDyn';
dynVideoFilename = join([dynvideoname p.run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
%lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
figure();
% hold on;
for t = 1:500:length(psi(1,1,:))-1
    plot1 = contourf(abs(dpsi1(:,:,t)),'levels',0.001, 'Linecolor', 'none');
    hold on;
    plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 2);
    hold off;
    set(gca,'YDir','normal');
    clim([0 0.1]);
    xlim([1 p.Nx]);
    ylim([1 p.Ny]);
    colorbar;
    colormap jet;
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    clear plot1 plot2;
    hold off;
end
hold off;
close(lat_dyn_video);

%% dpsi2 video with rolls
altframes = 2;
dynvideoname = 'gshrollSPVDyn';
dynVideoFilename = join([dynvideoname p.run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
%lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
figure();
% hold on;
for t = 1:500:length(psi(1,1,:))-1
    plot1 = contourf(abs(dpsi2(:,:,t)),'levels',0.001, 'Linecolor', 'none');
    hold on;
    plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 2);
    hold off;
    set(gca,'YDir','normal');
    clim([0 0.1]);
    xlim([1 p.Nx]);
    ylim([1 p.Ny]);
    colorbar;
    colormap jet;
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    clear plot1 plot2;
    hold off;
end
hold off;
close(lat_dyn_video);


%% dpsi3 video with rolls
altframes = 2;
dynvideoname = 'gshrollTPVDyn';
dynVideoFilename = join([dynvideoname p.run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
%lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
figure();
% hold on;
for t = 1:500:length(psi(1,1,:))-1
    plot1 = contourf(abs(dpsi3(:,:,t)),'levels',0.001, 'Linecolor', 'none');
    hold on;
    plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 2);
    hold off;
    set(gca,'YDir','normal');
    clim([0 0.1]);
    xlim([1 p.Nx]);
    ylim([1 p.Ny]);
    colorbar;
    colormap jet;
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    clear plot1 plot2;
    hold off;
end
hold off;
close(lat_dyn_video);

%% psi figure

% figure;
% imagesc(psimat(:,:,tmax)');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;
% title('\psi');

% figure;
% imagesc(ulat(:,:,6)');
% set(gca,'YDir','normal');
% colorbar;
% colormap jet;

%% rolls with zeta
figure;
contourf(psi(:,:,end-1),'levels',0.1, 'Linecolor', 'none');
clim([-1 1]);
hold on;
[X,Y] = meshgrid(1:p.Nx,1:p.Ny);
quiver(X,Y,v(:,:,end),u(:,:,end),2,'black');
hold off;
set(gca,'YDir','normal');
xlim([1 p.Nx]);
ylim([4 p.Ny]);
colorbar;
colormap jet;
title('\psi');

%% rolls without zeta
figure;
contourf(psi(:,:,end-1),'levels',0.1, 'Linecolor', 'none');
clim([-1 1]);
set(gca,'YDir','normal');
xlim([1 p.Nx]);
ylim([4 p.Ny]);
colorbar;
colormap jet;
title('\psi');

%%
figure;
contourf(dpsi1(:,:,end-1),'levels',0.005, 'Linecolor', 'none');
clim([-0.08 0.08]);
colormap jet;
title('\delta\psi_1');
colorbar;

%%
figure;
contourf(abs(dpsi1(:,:,end-1)),'levels',0.001, 'Linecolor', 'none');
hold on;
contourf(psi(:,:,end-1),'levels',1, 'Linecolor', 'black', ...
    'Facecolor', 'none', 'LineWidth', 2);
hold off;
set(gca,'YDir','normal');
clim([0 0.08]);
xlim([1 p.Nx]);
ylim([1 p.Ny]); 
colorbar;
colormap jet;
title('\psi and \delta\psi_1');

%%
figure;
contourf(omz(:,:,end-1),'levels',0.1, 'Linecolor', 'none');
clim([-1 1]);
colormap jet;
title('\Omega_z');
colorbar;

%%
figure;
contourf(abs(domz1(:,:,end-1)),'levels',0.005, 'Linecolor', 'none');
clim([0 0.1]);
colormap jet;
title('\delta\Omega_{z,1}');
colorbar;

%% 
figure;
plot(cumsum(lam1inst)./(1:length(lam1inst)), '-o');
xlabel('nmax/normtime');
ylabel('\lambda_{1,i}');

%% zeta contour figure

figure;
contourf(zeta(:,:,p.totimeu),'Linecolor','none');
set(gca,'YDir','normal');
colormap jet;
colorbar;
title('\zeta');


%% zeta contour figure with rolls
% timestep = tmax; % change this to desired time step (where spirals are seen)
% figure;
% p1 = contourf(abs(zetamat(:,:,timestep)),'levels',0.1, 'Linecolor', 'none');
% set(gca,'YDir','normal');
% xlim([1 Nx]);
% ylim([1 Ny]);
% colormap jet;
% colorbar;
% hold on;
% p2 = quiver(X,Y,2*vlat(:,:,timestep),2*ulat(:,:,timestep),2,'black');
% % clim([0 5]);
% p3 = contourf(psimat(:,:,timestep), 'levels', 1, 'Linecolor', 'black', ...
%     'Linewidth', 3,'Facecolor', 'none');
% hold off;

