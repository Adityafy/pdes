%% This script generates conditions for given parameters and save (for GSH)
% attempt to combine FDSI and PSETD run through a single script

clear all;
close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% Parameters

%% Control Parameters
con = initControlParams('csq',0);

%% SIMULATION PARAMETERS
sim.lam_0 = 2*pi; % critical roll wavelength (not a sim parameter but for now)
sim.Gamma = 40; % aspect ratio
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
sim.makeLiveFig = 0;
sim.saveTrDyn = 0;
sim.noSavingDynInLoop = 1;
sim.saveICforTS = 1;


%% ALL PARAMETERS:
p = paraStructGSH(con,sim,@spatialGridReal,@spatialGridSpectral,...
    @LinearOperatorGSHspectral);

csq = [0 0.2:0.2:1 2 4];
P = repmat(p, 1, length(csq));   % make a struct array

% Assign unique csq to each struct
for i = 1:length(csq)
    % P(i).con.csq = csq(i);
    % P(i).con.c = sqrt(csq(i));
    con = initControlParams('csq',csq(i));
    % optionally also update run name so files don’t overwrite
    P(i) = paraStructGSH(con,sim,@spatialGridReal,@spatialGridSpectral,...
    @LinearOperatorGSHspectral);
end

parfor i = 1:length(csq)
    fprintf('Running transients %s ...\n', P(i).dyn_run_name);

    [psiICpostTr, omzICpostTr, zetaICpostTr] = gshTimeIntgICpostTr(P(i), @gsh_IC_for_transients);

    fname = fullfile(sfa, [P(i).dyn_run_name 'ICpostTr.mat']);
    saveinparallel(fname, P(i), psiICpostTr, omzICpostTr, zetaICpostTr);
end


function saveinparallel(fname, p, psiICpostTr, omzICpostTr, zetaICpostTr)
    save(fname, 'p', 'psiICpostTr', 'omzICpostTr', 'zetaICpostTr', '-v7.3');
end


