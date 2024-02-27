
clear all;
close all;
addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address


%% Parameter struct

% key parameters
eps = 0.7;
sig = 2;
csq = 0.1;
gm = 50;
lam_0 = 2*pi;
dx = lam_0/8;
dy = dx;
dt = 0.001;
rolls = 5;
trtu = 20; % transient time units
totu = 100; % total time units
seed = 1;
normtime = 10;

p = paramgsh(eps,sig,csq,gm,dx,dy,dt,rolls,trtu,totu,normtime,seed);

%% Dynamics
tic
[psi, omz, zeta, ulat, vlat] = gshTimeIntgn(p,@gshtric);
toc

%% First Perturbation Vector
p.fpvictype = -1;
[dpsi1, domz1, dzeta1, fpv, fpvmag, lam1inst, lam1, res] = ...
    gshFpvExp(p,psi,omz,zeta,@rk2tsgsh1,@rk2tsgsh2,@gshfpvic);
toc

%%
% run('postprocessgsh.m');