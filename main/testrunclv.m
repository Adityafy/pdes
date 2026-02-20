%% CLV test run

clear all;
close all;

%%
addpath('../src/');

%%
sfa = '../../pdesDataDump/'; % saved files address

%% Load initial condtions from transients.
dynICAddress = join([sfa 'psG40eps0-7sig1csq0-1gm50tu100000s1ICpostTr']);
% load(dynICFileAddress);
load(dynICAddress,'p'); % load p to change some fields for ts calcs

%% Set up simulation parameters for tangent space calculations.
% we change the fields in simulation struct (sim) of p struct
dt_ts = 0.1;
tu_ts = 100;
p.sim.dt = dt_ts;
p.sim.tu = tu_ts;
p.sim.nmax = round(p.sim.tu/p.sim.dt);

%%%------saving type
%%% savingtype -> 0 means no saving big files at all
%%%            -> 1 means saving by parts but all time steps
%%%            -> 2 means saving only at defined time units ...
%%%                 ... (not all time steps, but still memory intensive)
%%%            -> 3 means save everything (ONLY FOR TESTING)
savingtype = 3;

bufInterv = 1;
allTUinterv = floor(p.sim.tu);
if p.sim.tu < 1
    allTUinterv = 10;
end
p.sim.savingtype = savingtype;
p.sim.interv = allTUinterv;
p.sim.bufInterv = bufInterv;
p.sim.makeLiveFig = 0; % 0/1 for no/yes figure on every time step
p.sim.progReportFactor = 10;

%% Set up etd for tangent space calculations (crucial if dt_ts is different)
p.etd = etdPreCalcs(p.L1,p.L2,p.sim.dt); % imperetive to update this

%% ts calculation parameters for reorthonormalization
tN = 1; % renormalization time units
nnorm = tN/p.sim.dt; % renormalization time step interval
nv = 2; % number of vectors to be calculated
ts = struct('tN',tN,'nnorm',nnorm,'nv',nv); % struct for tangent space parameters

% make the parameter struct from transients
p.ts = ts;

%%
fprintf(join(['Running tangent space ' p.dyn_run_name ' ...\n']));

%% Time Integration
tic;
% [psi,omz,zeta,dH,R,dHmag,laminst,lamgs] = gshTSTimeIntg(p,dynICAddress);
[psi,omz,zeta,dH,R,dHmag,laminst,lamgs] = gshTSTimeIntg(p,dynICAddress);
% [dHmag,laminst,lamgs] = gshTSTimeIntgBuffer(p,dynICAddress);
toc;

%% calculating cmatrix
% fprintf('\nCalcs for c matrix..\n');
% cmat = cmatrix(R,0,0);
% toc;
% clear R;

%%
fprintf(join(['Running clv calcs ' p.dyn_run_name ' ...\n']));

% [wf,cle] = clvGSH(p,psi,omz,zeta,dH,cmat);

[clv,cle,clvmag] = computeCLVs(p,dH,R,0,0);
toc;

%%
save(join([sfa p.dyn_run_name 'TS' num2str(p.sim.tu) 'nv' num2str(nv) 'clv']),'-v7.3');