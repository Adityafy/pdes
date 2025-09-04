%% Script for GSH 
% attempt to combine FDSI and PSETD run through a single script

clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% Parameters

%% Control Parameters
con = initControlParams('csq',0.1,'gm',50);

%% SIMULATION PARAMETERS
sim.lam_0 = 2*pi; % critical roll wavelength (not a sim parameter but for now)
sim.Gamma = 28; % aspect ratio
sim.dt = 0.1;
sim.tu = 100000; % time units
sim.nmax = round(sim.tu/sim.dt);
sim.seed = 1;
sim.frgamma = pi/2; %filtering radius gamma
sim.wlmult = sim.Gamma/2; % wavelength multiplier (determines Lx) (do not change)
sim.dynICtype = 1; % initial condition type for starting transients simulation
sim.runtype = 1; % 0 for finite difference + semi-Implicit
                % 1 for pseudospectral
sim.interv = floor(sim.tu); %intervals for saving dynamics

%%% switches
sim.progReportFactor = 10;
sim.makeLiveFig = 1;
sim.saveTrDyn = 0;
sim.noSavingDynInLoop = 1;
sim.saveICforTS = 1;


%% ALL PARAMETERS:
p = paraStructGSH(con,sim,@spatialGridReal,@spatialGridSpectral,...
                    @LinearOperatorGSHspectral);

fprintf(join(['Running transients ' p.dyn_run_name ' ...\n']));

%% Time Integration
tic;

if sim.saveTrDyn == 1
    [psitr,omztr,zetatr,utr,vtr] = gshTimeIntg(p,@gsh_IC_for_transients);

    psiICpostTr = psitr(:,:,end);
    omzICpostTr = omztr(:,:,end);
    zetaICpostTr = zetatr(:,:,end);
    uICpostTr = utr(:,:,end);
    vICpostTr = vtr(:,:,end);

    save(join([sfa p.dyn_run_name 'Transients']),'p','psitr','omztr','zetatr', ...
        'utr','vtr','psiICpostTr','omzICpostTr', 'zetaICpostTr', ...
        'uICpostTr','vICpostTr','-v7.3');
end

if sim.noSavingDynInLoop ==1
    [psiICpostTr, omzICpostTr, zetaICpostTr] = gshTimeIntgICpostTr(p, @gsh_IC_for_transients);
end

if sim.saveICforTS==1
    save(join([sfa p.dyn_run_name 'ICpostTr']),'p','psiICpostTr', ...
        'omzICpostTr', 'zetaICpostTr', 'uICpostTr','vICpostTr','-v7.3');
end


