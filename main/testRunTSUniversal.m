%% Script for GSH TS 
% attempt to combine FDSI and PSETD run through a single script

clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% Specify if this run is restarting from another tangent space file or after transients
isRestart = 1;

if isRestart ==1
    %  put tangent space run file here
    dynICAddress = join([sfa 'psG42eps0-7sig1csq0-1gm50time0-10s1TS10nv2']);
elseif isRestart ==0
    % put ICpostTr file address here
    dynICAddress = join([sfa 'psG42eps0-7sig1csq0-1gm50tu100000s1ICpostTr']);
end

% load parameter struct
load(dynICAddress,'p'); % load p to change some fields for ts calcs

% add a new restart field to p
p.Restart = isRestart;

%% Set up simulation parameters for tangent space calculations.

% change the fields in simulation struct (sim) of p struct
dt_ts = 0.1; % just there to be changed if necessary
tu_ts = 10; % MAKE SURE TO SPECIFY IT for this particular run
p.sim.dt = dt_ts;
p.sim.tu = tu_ts;
p.sim.nmax = round(p.sim.tu/p.sim.dt);

% new time start and end when restarted (ESSENTIAL IF RESTARTED)
if isRestart ==1
    % define new start time from end to new
    p.sim.tstart = p.sim.tend;
elseif isRestart ==0
    % starts from zero
    p.sim.tstart = 0;
end
p.sim.tend = p.sim.tstart + p.sim.tu;


%%%------saving type
%%% savingtype -> 0 means no saving big files at all
%%%            -> 1 means saving by parts but all time steps
%%%            -> 2 means saving only at defined time units ...
%%%                 ... (not all time steps, but still memory intensive)
savingtype = 2;

bufInterv = 1;
allTUinterv = floor(p.sim.tu)./2;
if p.sim.tu < 1
    allTUinterv = 10;
end

p.sim.savingtype = savingtype;
p.sim.interv = allTUinterv;
% p.sim.interv = allTUinterv/2;
p.sim.bufInterv = bufInterv;
p.sim.makeLiveFig = 1; % 0/1/2 for no/dynamics/cumulative lam1inst figure on every time step
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


%% update the run name
dyn_run_name = runName(p);
p.dyn_run_name = dyn_run_name;

%%
fprintf(join(['Running tangent space ' p.dyn_run_name ' ...\n']));

%% Time Integration
tic;
if isRestart ==1
    fprintf('This is a restart ...\n');
[psi,omz,zeta,dH,Rmat,dHmag,laminst,lamgs] = gshTSTimeIntgRestart(p,dynICAddress);
elseif isRestart ==0
[psi,omz,zeta,dH,dHmag,laminst,lamgs] = gshTSTimeIntg(p,dynICAddress);
% [dHmag,laminst,lamgs] = gshTSTimeIntgBuffer(p,dynICAddress);
end
toc;

%%
save(join([sfa p.dyn_run_name 'TS' num2str(p.sim.tu) 'nv' num2str(nv)]),'-v7.3');