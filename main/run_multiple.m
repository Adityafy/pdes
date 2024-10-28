
clear all;
close all;
addpath('../src/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

% lambda1 = [];

tic

csquare = [0.2:0.2:1 2 4];
% csquare = 0.2:0.2:0.4;
parfor i = 1:length(csquare)
    %% Parameter struct

    % key parameters
    eps = 0.7;
    sig = 1;
    csq = csquare(i);
    gm = 50;
    lam_0 = 2*pi; % critical roll wavelength
    dx = lam_0/8;
    dy = dx;
    dt_tr = 0.2;
    dt_fpv = 0.001;
    gamma = 15;
    trtu = 2000; % transient time units
    totu = 250; % total time units
    seed = 1;
    tN = 2; % renormalization time units

    p = paramgsh(eps,sig,csq,gm,lam_0,dx,dy,dt_tr,dt_fpv, ...
        gamma,trtu,totu,tN,seed);

    tic
    %% Transients
    [psitr, omztr, zetatr, utr, vtr] = gshTrSemiImpAtIntervals(p);

    % saving
    parsave1('dynonly', p, psitr, omztr, zetatr, utr, vtr);

    %% FPV calcs
    p.fpvictype = -1;
    [psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
        fpv, fpvmag, lam1inst, lam1, res] ...
        = gshFpvExpAtIntervals(p,psitr,omztr,zetatr);


    % saving
    parsave2('FPV', p, psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
        fpv, fpvmag, lam1inst, lam1, res);

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

function parsave1(name, p, psitr, omztr, zetatr, utr, vtr)
save(join([name p.run_name '.mat']),'psitr','omztr','zetatr','utr','vtr', '-v7.3');
end

function parsave2(name, p, psi, omz, zeta, u, v, dpsi1, domz1, dzeta1, ...
    fpv, fpvmag, lam1inst, lam1, res)
save(join([name p.run_name '.mat']), "p", "psi", "omz", "zeta", "u", "v", ...
    "dpsi1", "domz1", "dzeta1", ...
        "fpv", "fpvmag", "lam1inst", "lam1", "res", '-v7.3');
end