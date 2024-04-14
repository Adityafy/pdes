
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
rolls = 20;
trtu = 100; % transient time units
totu = 1000; % total time units
seed = 4;
normtime = 10;
intervals = 1000;

p = paramgsh(eps,sig,csq,gm,dx,dy,dt,rolls,trtu,totu,normtime,seed,intervals);

tic
%% Dynamics + FPV captured at intervals
p.fpvictype = -1;
[psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
                    fpv, fpvmag, lam1inst, lam1, res] = gshFpvExpAtIntervals(p);
toc

% %% Dynamics captured at intervals
% p.fpvictype = -1;
% tic
% [psi, omz, zeta, ulat, vlat] = gshTimeIntgnAtIntervals(p,10);
% toc
% 
% %% Dynamics
% tic
% [psi, omz, zeta, ulat, vlat] = gshTimeIntgn(p,@gshtric);
% toc
% 
% %% First Perturbation Vector
% p.fpvictype = -1;
% [dpsi1, domz1, dzeta1, fpv, fpvmag, lam1inst, lam1, res1] = ...
%     gshFpvExp(p,psi,omz,zeta,@rk2tsgsh1,@rk2tsgsh2,@gshfpvic);
% toc
% 
% %% Second Perturbation Vector
% p.fpvictype = -1;
% [dpsi2, domz2, dzeta2, spv, spvmag, lam2inst, lam2, res2] = ...
%     gshSpvExp(p,psi,omz,zeta,fpv,@rk2tsgsh1,@rk2tsgsh2,@gshfpvic);
% toc
% 
% %% Third Perturbation Vector
% p.fpvictype = -1;
% [dpsi3, domz3, dzeta3, tpv, tpvmag, lam3inst, lam3, res3] = ...
%         gshTpvExp(p,psi,omz,zeta,fpv,spv,@rk2tsgsh1,@rk2tsgsh2,@gshfpvic);
% toc

%% saving
save(join([p.run_name '.mat']), '-v7.3');

%%
% run('postprocessgsh.m');