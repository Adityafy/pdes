clear all;
close all;
addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address


%% Parameter struct

% key parameters
eps = 0.7;
sig = 1;
csq = 0.1;
gm = 50;
lam_0 = 2*pi; % critical roll wavelength
dx = lam_0/8;
dy = dx;
dt_tr = 0.2;
dt_fpv = 0.001;
gamma = 10;
trtu = 1000; % transient time units
totu = 200; % total time units
seed = 4;
tN = 2; % renormalization time units

p = paramgsh(eps,sig,csq,gm,lam_0,dx,dy,dt_tr,dt_fpv, ...
    gamma,trtu,totu,tN,seed);



tic
%% Transients
[psitr, omztr, zetatr, utr, vtr] = gshTrSemiImpAtIntervals(p);

save(join(['dynonly' p.run_name '.mat']), '-v7.3');

%% FPV calcs
p.fpvictype = -1;
[psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
    fpv, fpvmag, lam1inst, lam1, res] ...
        = gshFpvExpAtIntervals(p,psitr,omztr,zetatr);


%% saving
save(join([p.run_name '.mat']), '-v7.3');

%%
% run('postprocessgsh.m');