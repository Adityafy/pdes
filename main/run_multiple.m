
clear all;
close all;
addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

lambda1 = [];

tic

csquare = [0.2:0.2:1 2 4];
parfor i = 1:length(csquare)
    %% Parameter struct

    % key parameters
    eps = 0.7;
    sig = 2;
    csq = csquare(i);
    gm = 50;
    lam_0 = 2*pi;
    dx = lam_0/8;
    dy = dx;
    dt = 0.001;
    rolls = 1;
    trtu = 1; % transient time units
    totu = 1; % total time units
    seed = 3;
    normtime = 10;

    p = paramgsh(eps,sig,csq,gm,dx,dy,dt,rolls,trtu,totu,normtime,seed);
    [psi, omz, zeta, ulat, vlat] = gshTimeIntgn(p,@gshtric);
    toc

    %% First Perturbation Vector
    p.fpvictype = -1;
    [~, ~, ~, ~, ~, ~, lam1, ~] = ...
        gshFpvExp(p,psi,omz,zeta,@rk2tsgsh1,@rk2tsgsh2,@gshfpvic);
    toc
    
    lambda1 = [lambda1 lam1];
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
    % save(join([p.run_name '.mat']), '-v7.3');
end
%%
% run('postprocessgsh.m');