
clear all;
close all;
load('dynonly56x56eps0-7sig1csq0-1gm50tN2trt1000tot250s4.mat');
addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

tic
%% key parameters
eps = 0.7;
sig = 1;
csq = 0.2;
gm = 50;
lam_0 = 2*pi; % critical roll wavelength
dx = lam_0/8;
dy = dx;
dt_tr = 0.2;
dt_fpv = 0.001;
gamma = 15;
trtu = 2000; % transient time units
totu = 125; % total time units
seed = 4;
tN = 2; % renormalization time units

p = paramgsh(eps,sig,csq,gm,lam_0,dx,dy,dt_tr,dt_fpv, ...
    gamma,trtu,totu,tN,seed);


%% FPV

p.fpvictype = -1;
[psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
    fpv, fpvmag, lam1inst, lam1, res] ...
        = gshFpvExpAtIntervals(p,psitr,omztr,zetatr);

toc
%% saving
save(join([p.run_name '.mat']), '-v7.3');