%% Script for GSH 
% attempt to combine FDSI and PSETD run through a single script

clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% Parameters

% CONTROL PARAMETERS:
epsilon = 0.7;
sig = 1;
csq = 0.1;
c = sqrt(csq);
gm = 50;
% control parameters struct
con = struct('epsilon',epsilon,'sig',sig,'csq',csq,'c',c,'gm',gm); 

% SIMULATION PARAMETERS
lam_0 = 2*pi; % critical roll wavelength (not a sim parameter but for now)
dt = 0.1;
tu = 50; % time units
nmax = round(tu/dt);
seed = 1;
frgamma = pi/2; %filtering radius gamma
wlmult = 8; % wavelength multiplier (determines Lx)
dynICtype = 1; % initial condition type for starting transients simulation
runtype = 0; % 0 for finite difference + semi-Implicit
                % 1 for pseudospectral
interv = floor(tu); %intervals for saving dynamics
makeLiveFig = 1;
% simulation parameter struct
sim = struct('lam_0',lam_0,'dt',dt,'tu',tu,'seed',seed, 'nmax',nmax, ...
            'frgamma',frgamma,'wlmult',wlmult,'dynICtype',dynICtype, ...
            'runtype',runtype,'interv',interv,'makeLiveFig',makeLiveFig);

% ALL PARAMETERS:
p = paraStructGSH(con,sim,@spatialGridReal,@spatialGridSpectral,...
                    @LinearOperatorGSHspectral);


%% Intervals
% interv = floor(tu); %intervals
% TsIntervals = floor(tstu);
% if tstu < 1
%     TsIntervals = 10;
% end

fprintf(join(['Running transients ' p.dyn_run_name ' ...\n']));


%% For Finite-difference based GMRES (if used)
%%% making the finite difference matrix A (fdA) for Ax=b
% fdA = makeFDMMatrix(p);
restart = 20;
tol = 1e-6; % error tolerance
maxit = 200; % max # of iterations

%% Time Integration
tic;

[psitr,omztr,zetatr,utr,vtr] = gshTimeIntg(p,@gsh_IC_for_transients);

psiICpostTr = psitr(:,:,end);
omzICpostTr = omztr(:,:,end);
zetaICpostTr = zetatr(:,:,end);
uICpostTr = utr(:,:,end);
vICpostTr = vtr(:,:,end);


%% wavenumber figure
% figure;
% contourf(Kxn,Kyn,psihat,'LevelStep',0.0005,'EdgeColor','none');
% colorbar; colormap jet;
% axis square;
% set(gca,'TickLabelInterpreter','tex','FontSize',15);
% box("on");
% xlabel('$k_x$','Interpreter','latex','FontSize',30);
% ylabel('$k_y$','Interpreter','latex','FontSize',30);

%% saving transients
save(join([sfa p.dyn_run_name 'Transients']),'p','psitr','omztr','zetatr', ...
    'utr','vtr','psiICpostTr','omzICpostTr', 'zetaICpostTr', ...
    'uICpostTr','vICpostTr','-v7.3');
