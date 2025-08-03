%% Script for GSH TS 
% attempt to combine FDSI and PSETD run through a single script

clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% Load initial condtions from transients.
dynICAddress = join([sfa 'ps51N64eps0-7sig1csq0-1gm50tu2000s1Transients']);
% load(dynICFileAddress);
load(dynICAddress,'p');

%% Parameters for tangent space calculations.

% update the sim and etd field in p for tangent space calculations
dt_ts = 0.1;
tu_ts = 400;
p.sim.dt = dt_ts;
p.sim.tu = tu_ts;
p.sim.nmax = round(p.sim.tu/p.sim.dt);
interv = floor(p.sim.tu);
if p.sim.tu < 1
    interv = 10;
end
p.sim.interv = interv;
p.sim.makeLiveFig = 0; % 0/1 for no/yes figure on every time step
p.etd = etdPreCalcs(p.L1,p.L2,p.sim.dt);

% ts calculation parameters for reorthonormalization
tN = 1;
nnorm = tN/p.sim.dt;
nv = 3;
ts = struct('tN',tN,'nnorm',nnorm,'nv',nv);

% make the parameter struct from transients
p.ts = ts;

%%
fprintf(join(['Running tangent space ' p.dyn_run_name ' ...\n']));

%% Time Integration
tic;
[psi,omz,zeta,u,v,pertvecs,dHmag,laminst,lamgs] = gshTSTimeIntg(p,dynICAddress);

%%
save(join([sfa p.dyn_run_name 'TS' num2str(p.sim.tu) 'nv' num2str(nv)]),'-v7.3');