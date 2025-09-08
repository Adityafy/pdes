%% Script for GSH TS 
% attempt to combine FDSI and PSETD run through a single script

clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% Load initial condtions from transients.
dynICAddress = join([sfa 'psG30eps0-7sig1csq0-8gm50tu100000s1ICpostTr']);
% load(dynICFileAddress);
load(dynICAddress,'p'); % load p to change some fields for ts calcs

%% Set up simulation parameters for tangent space calculations.
% we change the fields in simulation struct (sim) of p struct
dt_ts = 0.1;
tu_ts = 900;
p.sim.dt = dt_ts;
p.sim.tu = tu_ts;
p.sim.nmax = round(p.sim.tu/p.sim.dt);

%%%------saving type
%%% savingtype -> 0 means no saving big files at all
%%%            -> 1 means saving by parts but all time steps
%%%            -> 2 means saving only at defined time units ...
%%%                 ... (not all time steps, but still memory intensive)
savingtype = 0;

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
tN = 2; % renormalization time units
nnorm = tN/p.sim.dt; % renormalization time step interval
nv = 1; % number of vectors to be calculated
ts = struct('tN',tN,'nnorm',nnorm,'nv',nv); % struct for tangent space parameters

% make the parameter struct from transients
p.ts = ts;

%%
fprintf(join(['Running tangent space ' p.dyn_run_name ' ...\n']));

%% Time Integration
tic;
[psi,omz,zeta,dH,dHmag,laminst,lamgs] = gshTSTimeIntg(p,dynICAddress);
% [dHmag,laminst,lamgs] = gshTSTimeIntgBuffer(p,dynICAddress);
toc;

%%
save(join([sfa p.dyn_run_name 'TS' num2str(p.sim.tu) 'nv' num2str(nv)]),'-v7.3');