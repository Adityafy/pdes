
addpath('../../../pdes/src/');
addpath('../../../pdes/src/othercolor/');

%% Transient dynamics video

altframes = 2;

    dynvideoname = 'gshTrDyn';
    dynVideoFilename = join([dynvideoname p.run_name]);
    lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
    lat_dyn_video.FrameRate = 30;
    open(lat_dyn_video);
    [X,Y] = meshgrid(1:p.Nx,1:p.Ny);
    figure('units','pixels','position',[0 0 1440 1080]);
    
    % hold on;
    for t = 1:1:length(psitr(1,1,:))-1
        % imagesc(psilat(:,:,t));
        plot1 = contourf(psitr(:,:,t),'levels',0.1, 'Linecolor', 'none');
        set(gca,'YDir','normal');
        hold on;
        % imagesc(zetalat(:,:,t));
        colorbar;
        colormap jet
        clim([-0.8 0.8]);
        plot2 = quiver(X,Y,vtr(:,:,t),utr(:,:,t),2,'black');
        xlim([1 p.Nx]);
        ylim([1 p.Ny]);
        axis square;
        frame = getframe(gcf);
        writeVideo(lat_dyn_video,frame);
        clear plot1 plot2;
        hold off;
    end
    hold off;
    close(lat_dyn_video);

%% Dynamics video (when no p present)

dynvideoname = 'gshDyn';
dynVideoFilename = join([dynvideoname p.dyn_run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
x = p.rmesh.x; y = p.rmesh.y; Nx = p.rmesh.Nx; Ny = p.rmesh.Ny;
figure;
% figure('units','pixels','position',[0 0 1440 1080]);

% hold on;
for t = 1:1:length(psi(1,1,:))
    % imagesc(psilat(:,:,t));
    plot1 = contourf(X,Y,psi(:,:,t),'levels',0.1, 'Linecolor', 'none');
    set(gca,'YDir','normal');
    hold on;
    % imagesc(zetalat(:,:,t));
    colorbar;
    colormap jet
    clim([-0.8 0.8]);
    % plot2 = quiver(X,Y,u(:,:,t),v(:,:,t),2,'black');
    % xlim([1 Nx]);
    % ylim([1 Ny]);
    axis square;
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    clear plot1 plot2;
    hold off;
end
hold off;
close(lat_dyn_video);

%% dynamics video with mean flow vector arrows

altframes = 2;

    dynvideoname = 'gshDyn';
    dynVideoFilename = join([dynvideoname p.run_name]);
    lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
    lat_dyn_video.Quality = 100;
    %lat_dyn_video.FrameRate = 30;
    open(lat_dyn_video);
    [X,Y] = meshgrid(1:p.Nx,1:p.Ny);
    figure();
    
    % hold on;
    for t = 1:1:length(psi(1,1,:))-1
        % imagesc(psilat(:,:,t));
        plot1 = contourf(psi(:,:,t),'levels',0.1, 'Linecolor', 'none');
        % plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        % 'Facecolor', 'none', 'LineWidth', 1.5);
        set(gca,'YDir','normal');
        hold on;
        % imagesc(zetalat(:,:,t));
        % colorbar;
        colormap jet;
        
        clim([-0.8 0.8]);
        colorbar(gca,'Ticks',[-0.8 0.8]);
        plot2 = quiver(X,Y,v(:,:,t),u(:,:,t),3,'black');
        xlim([1 p.Nx]);
        ylim([1 p.Ny]);
        axis square;
        frame = getframe(gcf);
        writeVideo(lat_dyn_video,frame);
        clear plot1 plot2;
        hold off;
    end
    hold off;
    close(lat_dyn_video);

    %% psi only video


    dynvideoname = 'psi';
    dynVideoFilename = join([dynvideoname p.run_name]);
    lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
    lat_dyn_video.Quality = 100;
    %lat_dyn_video.FrameRate = 30;
    open(lat_dyn_video);
    figure();
    
    % hold on;
    for t = 1:1:length(psi(1,1,:))-1
        plot1 = contourf(psi(:,:,t),'levels',0.1, 'Linecolor', 'none');
        set(gca,'YDir','normal');
        % colorbar;
        colormap jet;
        clim([-0.8 0.8]);
        colorbar(gca,'Ticks',[-0.8 0.8]);
        axis square;
        frame = getframe(gcf);
        writeVideo(lat_dyn_video,frame);
        clear plot1;
    end
    hold off;
    close(lat_dyn_video);

%% video : mean flow magnitude with rolls

altframes = 2;

    dynvideoname = 'mfs';
    dynVideoFilename = join([dynvideoname p.run_name]);
    lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
    %lat_dyn_video.FrameRate = 30;
    open(lat_dyn_video);
    [X,Y] = meshgrid(1:p.Nx,1:p.Ny);
    mfm = zeros(p.Nx,p.Ny,p.totimeu);
    for t = 1:1:length(psi(1,1,:))-1
        mfm(:,:,t) = sqrt(u(:,:,t).^2+v(:,:,t).^2);
    end

    figure();
    
    % hold on;
    for t = 1:1:length(psi(1,1,:))
        % mfm = sqrt(u(:,:,t).^2+v(:,:,t).^2);
        plot1 = contourf(mfm(:,:,t),'levels',0.05,'Linecolor','none');
        hold on;
        plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
                'Facecolor', 'none', 'LineWidth', 1.5);
        plot3 = quiver(X,Y,v(:,:,t),u(:,:,t),2,'black');
        
        colorbar;
        clim([0 2]);
        colormap jet
        xlim([1 p.Nx]);
        ylim([1 p.Ny]);
        axis square;

        frame = getframe(gcf);
        writeVideo(lat_dyn_video,frame);

        clear plot1 plot2;
        hold off;
    end
    hold off;
    close(lat_dyn_video);
%% plot : mean flow magnitude with rolls
[X,Y] = meshgrid(1:p.Nx,1:p.Ny);
mfm = zeros(p.Nx,p.Ny,p.totimeu);
for t = 1:1:length(psi(1,1,:))-1
    mfm(:,:,t) = sqrt(u(:,:,t).^2+v(:,:,t).^2);
end

figure;
attime = 955;
plot1 = contourf(mfm(:,:,attime),'levels',0.05,'Linecolor','none');
hold on;
plot2 = contourf(psi(:,:,attime),'levels',1, 'Linecolor', 'black', ...
    'Facecolor', 'none', 'LineWidth', 1.5);
plot3 = quiver(X,Y,v(:,:,attime),u(:,:,attime),2,'black');
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
set(gca,'YTickLabel',[],'XTickLabel',[]);
colorbar;
clim([0 2]);
colormap jet
xlim([1 p.Nx]);
ylim([1 p.Ny]);
axis square;
% subtitle '(a)';
hold off;

figure;
plot4 = contourf(abs(dpsi1(:,:,attime)),'levels',0.001, 'Linecolor', 'none');
hold on;
plot5 = contourf(psi(:,:,attime),'levels',1, 'Linecolor', 'black', ...
    'Facecolor', 'none', 'LineWidth', 1.5);

set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
set(gca,'YTickLabel',[],'XTickLabel',[]);
clim([0 0.05]);
xlim([1 p.Nx]);
ylim([1 p.Ny]);
axis square;
colorbar(gca,'Ticks',[0 0.05]);
colormap parula;
% colormap(colmap);
% subtitle '(b)';
hold off;

%% NEW dpsi1 video with rolls 
% N = p.rmesh.N;
% Nx = p.rmesh.Nx;
% Ny = p.rmesh.Ny;
dH1 = squeeze(pertvecs(:,1,:));
for n = 1:length(pertvecs(1,1,:))
    dpsi1(:,:,n) = reshape(pertvecs(1:p.rmesh.N,1,n),p.rmesh.Nx,p.rmesh.Ny)';
    domz1(:,:,n) = reshape(pertvecs(p.rmesh.N+1:2*p.rmesh.N,1,n),p.rmesh.Nx,p.rmesh.Ny)';
end

%%
for n = 1:length(dH1(1,:))
    dpsi1(:,:,n) = reshape(dH1(1:p.rmesh.N,n),p.rmesh.Nx,p.rmesh.Ny)';
    domz1(:,:,n) = reshape(dH1(p.rmesh.N+1:2*p.rmesh.N,n),p.rmesh.Nx,p.rmesh.Ny)';
end

%%
dynvideoname = 'dpsi1';
dynVideoFilename = join([dynvideoname p.dyn_run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
lat_dyn_video.Quality = 100;
%lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
figure;

% hold on;
for t = 1:1:length(psi(1,1,:))
    plot1 = contourf(abs(dpsi1(:,:,t)),'levels',0.001, 'Linecolor', 'none');
    hold on;
    plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 2);
    hold off;
    set(gca,'YDir','normal');
    set(gca,'YTickLabel',[],'XTickLabel',[]);
    clim([0 0.1]);
    xlim([1 p.rmesh.Nx]);
    ylim([1 p.rmesh.Ny]);
    axis square;
    colorbar;
    colormap jet;
    % colorbar(gca,'Ticks',[-0.005 0.05],'TickLabels',[0 0.05]);
    colorbar(gca,'Ticks',[0 0.1]);
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    clear plot1 plot2;
    hold off;
end
hold off;
close(lat_dyn_video);

%% NEW dpsi1,2,3 video with rolls 
N = p.rmesh.N;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;
% for n = 1:length(pertvecs(1,1,:))
%     dpsi1(:,:,n) = reshape(pertvecs(1:N,1,n),Nx,Ny)';
%     dpsi2(:,:,n) = reshape(pertvecs(1:N,2,n),Nx,Ny)';
%     dpsi3(:,:,n) = reshape(pertvecs(1:N,3,n),Nx,Ny)';
% end
dynvideoname = 'dH123';
dynVideoFilename = join([dynvideoname p.dyn_run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
lat_dyn_video.Quality = 100;
%lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
figure('units','pixels','position',[0 0 1440 1080]);
% hold on;
tlayout = tiledlayout(1, nv, 'Padding', 'compact', 'TileSpacing', 'compact'); % 20-ish px spacing
for t = 1:1:length(psi(1,1,:))
    for k = 1:p.ts.nv
        % subplot(1,p.ts.nv,k);
        nexttile(k);
        plot1 = contourf(abs(reshape(pertvecs(1:N,k,t),Nx,Ny)'),'levels',0.001, 'Linecolor', 'none');
        hold on;
        plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
            'Facecolor', 'none', 'LineWidth', 2);
        hold off;
        set(gca,'YDir','normal');
        set(gca,'YTickLabel',[],'XTickLabel',[]);
        title(join(['$\|\delta \psi^{(' num2str(k) ')} \|$']),'Interpreter','latex');
        clim([0 0.1]);
        xlim([1 p.rmesh.Nx]);
        ylim([1 p.rmesh.Ny]);
        axis square;
        colorbar;
        colormap jet;
        % colorbar(gca,'Ticks',[-0.005 0.05],'TickLabels',[0 0.05]);
        colorbar(gca,'Ticks',[0 0.1]);
        hold off;
    end
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    % clear plot1 plot2;
    tlayout = tiledlayout(1, nv, 'Padding', 'compact', 'TileSpacing', 'compact'); % 20-ish px spacing
end
hold off;
close(lat_dyn_video);

%%
N = p.rmesh.N;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;

% Video setup
dynvideoname = 'dH123';
dynVideoFilename = join([dynvideoname p.dyn_run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
lat_dyn_video.Quality = 100;
lat_dyn_video.FrameRate = 10;
open(lat_dyn_video);

nv = p.ts.nv;

fig = figure('units','pixels','position',[0 0 1440 1080]);
% fig = figure;
tlayout = tiledlayout(1, nv, 'Padding', 'compact', 'TileSpacing', 'compact'); % 20-ish px spacing

% Loop over time steps
for t = 1:size(psi, 3)
    for k = 1:p.ts.nv
        nexttile(k);
        contourf(abs(reshape(pertvecs(1:N,k,t),Nx,Ny)'), ...
            'levels',0.001, 'Linecolor', 'none');
        hold on;
        contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
            'Facecolor', 'none', 'LineWidth', 2);
        set(gca,'YDir','normal');
        set(gca,'YTickLabel',[],'XTickLabel',[]);
        clim([0 0.1]);
        xlim([1 p.rmesh.Nx]);
        ylim([1 p.rmesh.Ny]);
        axis square;
        colorbar;
        colormap jet;
        % colorbar(gca,'Ticks',[-0.005 0.05],'TickLabels',[0 0.05]);
        colorbar(gca,'Ticks',[0 0.1]);
        hold off;
    end

    % Export graphics to image and convert to frame
    exportgraphics(fig, 'frame.png');
    frameimg = imread('frame.png');
    [h, w, ~] = size(frameimg);
    h_even = 2 * floor(h / 2);
    w_even = 2 * floor(w / 2);
    frameimg = frameimg(1:h_even, 1:w_even, :);
    writeVideo(lat_dyn_video, im2frame(frameimg));

    % Optional: clear tiles to save memory
    clf(fig);
    tlayout = tiledlayout(1, p.ts.nv, 'Padding', 'tight', 'TileSpacing', 'none');
end

close(lat_dyn_video);
close(fig);


%% dpsi1 video with rolls
addpath('~/Documents/pdes/src/othercolor/');
colmap = othercolor('YlOrRd9');
altframes = 2;
dynvideoname = 'FPVorange';
dynVideoFilename = join([dynvideoname p.run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
lat_dyn_video.Quality = 100;
%lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
figure;

% hold on;
for t = 1:1:length(psi(1,1,:))-1
    plot1 = contourf(abs(dpsi1(:,:,t)),'levels',0.001, 'Linecolor', 'none');
    hold on;
    plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 2);
    hold off;
    set(gca,'YDir','normal');
    set(gca,'YTickLabel',[],'XTickLabel',[]);
    clim([0 0.05]);
    xlim([1 p.Nx]);
    ylim([1 p.Ny]);
    axis square;
    % colorbar;
    % colormap jet;
    % colorbar(gca,'Ticks',[-0.005 0.05],'TickLabels',[0 0.05]);
    colorbar(gca,'Ticks',[0 0.05]);
    colormap(colmap);
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    clear plot1 plot2;
    hold off;
end
hold off;
close(lat_dyn_video);

%% video of subplot = mean flow strength + first perturbation vector
addpath('~/Documents/pdes/src/othercolor/');
colmap = othercolor('YlOrRd9');

[X,Y] = meshgrid(1:p.Nx,1:p.Ny);
mfm = zeros(p.Nx,p.Ny,p.totimeu);
for t = 1:1:length(psi(1,1,:))
    mfm(:,:,t) = sqrt(u(:,:,t).^2+v(:,:,t).^2);
end

specificvidfilename = 'mfsfpv_col1_';
vidfilename = join([specificvidfilename p.run_name]);
video = VideoWriter(vidfilename, 'MPEG-4');
video.FrameRate = 10;
video.Quality = 100;
open(video);

vidfig = figure;

% set(gcf, 'Position',  [100, 100, 500, 400]);
for t = 1:1:length(psi(1,1,:))
    subplot(1,2,1);
    plot1 = contourf(mfm(:,:,t),'levels',0.05,'Linecolor','none');
    hold on;
    plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 1);
    plot3 = quiver(X,Y,v(:,:,t),u(:,:,t),2,'black');
    colorbar(gca,'Ticks',[0 2]);
    clim([0 2]);
    colormap(jet);
    xlim([1 p.Nx]);
    ylim([1 p.Ny]);
    set(gca,'YTickLabel',[],'XTickLabel',[]);
    axis square;
    hold off

    subplot(1,2,2);
    plot4 = contourf(abs(dpsi1(:,:,t)),'levels',0.001, 'Linecolor', 'none');
    hold on;
    plot5 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 1);

    set(gca,'YDir','normal');
    set(gca,'YTickLabel',[],'XTickLabel',[]);
    clim([0 0.05]);
    xlim([1 p.Nx]);
    ylim([1 p.Ny]);
    axis square;
    colorbar(gca,'Ticks',[0 0.05]);
    % colormap jet;
    colormap(colmap);
    hold off;
    
    frame = getframe(gcf);
    writeVideo(video,frame);
    clear plot1 plot2 plot3 plot4 plot5;
end
close(video);

%% dpsi1 video with rolls
altframes = 2;
dynvideoname = 'gshrollFPVDynfr10_';
dynVideoFilename = join([dynvideoname p.run_name]);
lat_dyn_video = VideoWriter(dynVideoFilename, 'MPEG-4');
lat_dyn_video.FrameRate = 10;
open(lat_dyn_video);
figure();
% hold on;
for t = 1:1:length(psi(1,1,:))-1
    plot1 = contourf(abs(dpsi1(:,:,t)),'levels',0.001, 'Linecolor', 'none');
    hold on;
    plot2 = contourf(psi(:,:,t),'levels',1, 'Linecolor', 'black', ...
        'Facecolor', 'none', 'LineWidth', 2);
    hold off;
    set(gca,'YDir','normal');
    clim([0 0.07]);
    xlim([1 p.Nx]);
    ylim([1 p.Ny]);
    axis square;
    colorbar;
    colormap jet;
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    clear plot1 plot2;
    hold off;
end
hold off;
close(lat_dyn_video);

%% dpsi1 histogram
lat_dyn_video = VideoWriter('dpsi1histogram', 'MPEG-4');
lat_dyn_video.FrameRate = 30;
open(lat_dyn_video);
figure();
% hold on;
for t = 1:1:length(psi(1,1,:))
    plot1 = histogram(latToVec(dpsi1(:,:,t)));
    axis square;
    xlim([-0.1 0.1]);
    ylim([0 1000]);
    frame = getframe(gcf);
    writeVideo(lat_dyn_video,frame);
    clear plot1;
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
    axis square;
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
    axis square;
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

%% psi rolls with zeta
figure;
contourf(psi(:,:,end),'levels',0.1, 'Linecolor', 'none');
clim([-1 1]);
hold on;
[X,Y] = meshgrid(1:p.Nx,1:p.Ny);
quiver(X,Y,v(:,:,end),u(:,:,end),2,'black');
hold off;
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
set(gca,'YTickLabel',[],'XTickLabel',[]);
xlim([1 p.Nx]);
ylim([4 p.Ny]);
axis square;
% axis off;
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
colorbar;
colormap jet;
title('\psi');

%% psi rolls with zeta and mean flow magnitude
figure1 = figure;
axes1 = axes('Parent',figure1);
hold on;
mfm = sqrt(u(:,:,end).^2+v(:,:,end).^2); % mean flow magnitude
contourf(mfm,'levels',0.05,'Linecolor','none');
contourf(psi(:,:,end),'levels',1, 'Linecolor', 'black', ...
    'Facecolor', 'none', 'LineWidth', 2);
[X,Y] = meshgrid(1:p.Nx,1:p.Ny);
quiver(X,Y,v(:,:,end),u(:,:,end),2,'black');

hold off;
set(gca,'YDir','normal');
set(axes1,'CLim',[0 2],'FontSize',15,'XTick',zeros(1,0),'YTick',zeros(1,0));
colormap(othercolor('Blues4'));
colorbar(axes1,'Ticks',[0 2]);
xlim([1 p.Nx]);
ylim([4 p.Ny]);
axis square;
title('$\left|\vec{U}\right|$','Interpreter','latex');

%% (NEW ps) psi rolls with zeta and mean flow magnitude
figure1 = figure;
X = p.rmesh.X;
Y = p.rmesh.Y;
totaltimesteps = size(psi,3);
fractn_howfar_from_end = 0.58;
at_n = totaltimesteps - round(fractn_howfar_from_end * totaltimesteps);
axes1 = axes('Parent',figure1);
hold on;
mfm = sqrt(u(:,:,at_n).^2+v(:,:,at_n).^2); % mean flow magnitude
contourf(X,Y,mfm,'levels',0.05,'Linecolor','none');
contourf(X,Y,psi(:,:,at_n),'levels',1, 'Linecolor', 'black', ...
    'Facecolor', 'none', 'LineWidth', 1.5);

quiver(X,Y,u(:,:,at_n),v(:,:,at_n),2,'black');

hold off;
set(gca,'YDir','normal');
set(axes1,'CLim',[0 3],'FontSize',15,'XTick',zeros(1,0),'YTick',zeros(1,0));
% colormap(othercolor('Blues4'));
colormap jet;
colorbar(axes1,'Ticks',[0 3]);
xlim([X(1,1) X(end,end)]);
ylim([Y(1,1) Y(end,end)]);
axis square;
title('$\left|\vec{U}\right|$','Interpreter','latex');

%% doublepsi
figure;
hold on;
% mfm = sqrt(u(:,:,end).^2+v(:,:,end).^2); % mean flow magnitude
contourf(psi(:,:,end),'levels',0.05,'Linecolor','none');
contourf(psi(:,:,end),'levels',1, 'Linecolor', 'black', ...
    'Facecolor', 'none', 'LineWidth', 1);
% [X,Y] = meshgrid(1:p.Nx,1:p.Ny);
% quiver(X,Y,v(:,:,end),u(:,:,end),2,'black');

hold off;
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
xlim([1 p.Nx]);
ylim([4 p.Ny]);
axis square;
colorbar;
clim([-1 1]);
colormap jet;
title('\psi');

%% psi rolls without zeta
figure;
contourf(psi(:,:,955),'levels',0.1, 'Linecolor', 'none');
clim([-1 1]);
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
set(gca,'YTickLabel',[],'XTickLabel',[]);
xlim([1 p.Nx]);
ylim([1 p.Ny]);
axis square;
colorbar;
colorbar(gca,'Ticks',[-1 1]);
colormap jet;
% title('\psi');

%% only dpsi1
figure;
contourf(dpsi1(:,:,end),'levels',0.005, 'Linecolor', 'none');
clim([-0.08 0.08]);
axis square;
colormap jet;
title('\delta\psi_1');
colorbar;

%% dpsi1 with rolls
figure;
contourf(abs(dpsi1(:,:,end)),'levels',0.001, 'Linecolor', 'none');
hold on;
contourf(psi(:,:,end),'levels',1, 'Linecolor', 'black', ...
    'Facecolor', 'none', 'LineWidth', 2);
hold off;
set(gca,'YDir','normal');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
clim([0 0.05]);
xlim([1 p.Nx]);
ylim([1 p.Ny]); 
axis square;
colorbar;
colormap jet;
title('\psi and \delta\psi_1');

%% omz
figure;
contourf(omz(:,:,end),'levels',0.1, 'Linecolor', 'none');
clim([-1 1]);
axis square;
colormap jet;
title('\Omega_z');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
colorbar;

%% domz1
figure;
contourf(abs(domz1(:,:,end)),'levels',0.005, 'Linecolor', 'none');
clim([0 0.1]);
axis square;
colormap jet;
title('\delta\Omega_{z,1}');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
colorbar;

%% lambda_1 cummulative average (lam1ca)
lamcafig = figure;
axlamfig = axes('Parent',lamcafig);
hold(axlamfig,'on');

hold on;
time = linspace(1,p.sim.tu,floor(p.sim.tu/p.ts.tN));
for k = 1:p.ts.nv
    lamca = cumsum(laminst(k,:))./(1:length(laminst(k,:)));
    plot(time,lamca,'MarkerSize',10,'Marker','.','LineWidth',1,'DisplayName',join(['k = ' num2str(k)]));
end
hold off;
% plot(cumsum(lam1inst)./(1:length(lam1inst)), '-o');
xlim([1 time(end)]);
ylim([0 1]);
% yline(0,'Parent',axlam1fig);
legend;

box(axlamfig,'on');
axis(axlamfig,'square');
hold(axlamfig,'off');
set(axlamfig,'FontSize',15,'LineWidth',1,'XMinorTick','on','YMinorTick','on');
ylabel('$\left<\lambda_{k}(t)\right>$','FontSize',30,'Interpreter','latex','Rotation',0);
xlabel('$t$','FontSize',30,'Interpreter','latex');

%% Lyapunov Spectrum
nvvec = (1:p.ts.nv)';
myplot(nvvec,lamgs,'$k$','$\lambda_k$');


%%
dlamb = kydimension(lamgs,length(lamgs));

%%
dlamfit = kydimension(lamfit,length(lamfit));

%% Fit: 'lam1 cummulative average fit'.
time = linspace(1,p.sim.tu,floor(p.sim.tu/p.ts.tN));
lam1ca = cumsum(laminst(1,:))./(1:length(laminst(1,:)));
fprintf('lambda_1 : %.4f\n\n',lamgs(1));
[xData, yData] = prepareCurveData( time, lam1ca );

fiteqn = 'a+b/x';
initExcludeXDataPoints = 50;

% Set up fittype and options.
ft = fittype( fiteqn, 'independent', 'x', 'dependent', 'y' );
excludedPoints = xData < initExcludeXDataPoints;
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.1 -1];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fprintf('Fitresult:\n');
disp(fitresult);
    
% Plot fit with data.
figure;
hold on;
h = plot( fitresult,'b-', xData, yData,'k.', excludedPoints, 'r*');
ylim([-max(abs(yData)) max(abs(yData))]);
xlim([xData(1) xData(end)]);
xaxis = yline(0);
set(h,'LineWidth',1,'MarkerSize',6);
fitname = join(['(' num2str(fitresult.a) ') + (' num2str(fitresult.b) ')/t']);
excludeLegend = join(['Initial points excluded: ' num2str(initExcludeXDataPoints)]);
legend( h, '$\left<\lambda_{1}(t)\right>_t/t$', excludeLegend, ...
    fitname, 'Location', 'SouthEast', 'Interpreter', 'latex' );
% Label axes
set(gca,'TickLabelInterpreter','tex','FontSize',15,'LineWidth',1,'XMinorTick','on','YMinorTick','on');
ylabel('$\left<\lambda_{1}(t)\right>_t/t$','FontSize',30,'Interpreter','latex','Rotation',90);
% ylabel('$(\sum_{i=0}^t\lambda_{1,\,i})/t$','FontSize',30,'Interpreter','latex','Rotation',0);
xlabel('$t$','FontSize',30,'Interpreter','latex');
box on; axis square;

%% Mean flow: maximum and threshold
% addpath('../../../pdes/src/');
maxmfst = zeros(1,p.totimeu);
mfm = zeros(p.Nx,p.Ny,p.totimeu);
mfthreshold = 1;
fpvthreshold = 0.01;
highmf = zeros(p.Nx,p.Ny,p.totimeu);
highmfratio = zeros(1,p.totimeu);
for t = 1:1:length(u(1,1,:))
    % mean flow magnitude
    mfm(:,:,t) = sqrt(u(:,:,t).^2+v(:,:,t).^2);
    maxmfst(t) = max(max(mfm(:,:,t)));
    
    highmf(:,:,t) = mfm(:,:,t) > mfthreshold;
    highdpsi1(:,:,t) = abs(dpsi1(:,:,t)) > fpvthreshold;
    corrmat_hfm_hfpv(:,:,t) = highmf(:,:,t).*highdpsi1(:,:,t);
    corrsumt(t) = sum(sum(corrmat_hfm_hfpv(:,:,t)))./(p.Nx*p.Ny);
    highmfratio(t) = sum(sum(highmf(:,:,t)))/(p.Nx*p.Ny);
    mfsum(t) = sum(sum(mfm(:,:,t)))/(p.Nx*p.Ny);
end
maxmfs = mean(maxmfst);
highmffracavg = mean(highmfratio);
highmffracstd = std(highmfratio);
mfmavg = sum(mfsum)/p.totimeu;
corravg = sum(corrsumt)/length(corrsumt);
corrsumstd = std(corrsumt);

diary meanflowresults;
diary on;
fprintf('\n--------- ');
fprintf(datestr(datetime('today')));
fprintf(' ----------\n');
fprintf(p.run_name);
fprintf('\n');
fprintf('\n Maximum mean flow magnitude: %.4f \n', maxmfs);
fprintf(['Time averaged percentage' ...
    ' of high mean flow: %.4f, (threshold: %g) \n'], highmffracavg,mfthreshold);
fprintf('Mean flow average : %.4f\n', mfmavg);
fprintf('-------------------------\n');
diary off;

%
save('corr_hfm_hfpv.mat','mfm','highdpsi1','highmf','corrmat_hfm_hfpv','highmfratio',...
'corrsumt','corravg','corrsumstd','maxmfs',"highmffracavg","mfmavg","highmffracstd");

% 
maxmfs_sweep = [maxmfs_sweep maxmfs];
highmffracavg_sweep = [highmffracavg_sweep highmffracavg];
highmffracstd_sweep = [highmffracstd_sweep highmffracstd];
mfmavg_sweep = [mfmavg_sweep mfmavg];
corravg_sweep = [corravg_sweep corravg];
corrsumstd_sweep = [corrsumstd_sweep corrsumstd];

clearvars -except maxmfs_sweep highmffracavg_sweep mfmavg_sweep ...
    corravg_sweep corrsumstd_sweep highmffracstd_sweep;

%% DO NOT RUN
% maxmfs_sweep = [];
% highmffracavg_sweep = [];
% mfmavg_sweep = [];
% corravg_sweep = [];
% corrsumstd_sweep = [];
% highmffracstd_sweep = [];

%% High mean flow and FPV
figure;
imagesc(mfm(:,:,end));
figure;
imagesc(dpsi1(:,:,end));

%% FPV magnitude
figure; plot(dHmag(1,:),'-o',Marker='.');
xlim([1 length(dHmag(1,:))]);
axis square;
xlabel('nmax');
ylabel('$||\delta \vec{m}||$','Interpreter','latex');
set(gca,'TickLabelInterpreter','tex','FontSize',15);

%% FPV magnitude
t1 = p.sim.tu - round(0.9*p.sim.tu) + 1;
t2 = p.sim.tu;
figure; plot(t1:t2,dHmag(1,t1:t2),'-o',Marker='.');
axis square;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
xlabel('$t$','Interpreter','latex','FontSize',30);
ylabel('$\|\delta \vec{H}^{(1)}\|$','Interpreter','latex', ...
    'Rotation',0,'FontSize',30);


%% zeta contour figure

figure;
contourf(zeta(:,:,p.sim.tu),'Linecolor','none');
set(gca,'YDir','normal');
axis square;
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
% axis square;
% colormap jet;
% colorbar;
% hold on;
% p2 = quiver(X,Y,2*vlat(:,:,timestep),2*ulat(:,:,timestep),2,'black');
% % clim([0 5]);
% p3 = contourf(psimat(:,:,timestep), 'levels', 1, 'Linecolor', 'black', ...
%     'Linewidth', 3,'Facecolor', 'none');
% hold off;
%%
Dlfit = kydimension(X,256);
%%

function D = kydimension(lamgs,M)
%
% D_lambda = kydimension(lamgs,M)
%
% Calculates the Kaplan-Yorke Dimension
%
%  lamgs    = vector of lyapunov exponents
%  M        = system size

jfd = 0; % jfd here is the j in Kaplan-Yorke Formula
for j = 1:M-1
    if sum(lamgs(1:j))>0 && sum(lamgs(1:j+1))<0
        jfd = j;
    end
end
D = jfd + sum(lamgs(1:jfd))/abs(lamgs(jfd+1));
if D == 0
    D = -1;
end

end
