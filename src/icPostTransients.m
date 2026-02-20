function icPostTransients(epsilon,sig,csq,gm,Gamma,transientTime)

%% IC post Transients, gamma = 38

% clear all;
% close all;

addpath('../src/');
sfa = '../../pdesDataDump/'; % saved files address

%% Parameters

%% Control Parameters
con = initControlParams('epsilon',epsilon,'sig',sig,'csq',csq,'gm',gm);
tu = transientTime;
%% SIMULATION PARAMETERS
sim.lam_0 = 2*pi; % critical roll wavelength (not a sim parameter but for now)
sim.Gamma = Gamma; % aspect ratio
sim.dt = 0.1;
sim.tu = tu; % time units
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
sim.saveTrDyn = 0; % DO NOT TURN THIS ON unless absolutely necessary
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
        'psiICpostTr','omzICpostTr', 'zetaICpostTr', '-v7.3');
end

if sim.noSavingDynInLoop ==1
    [psiICpostTr, omzICpostTr, zetaICpostTr] = gshTimeIntgICpostTr(p, @gsh_IC_for_transients);
end

if sim.saveICforTS==1
    save(join([sfa p.dyn_run_name 'ICpostTr']),'p','psiICpostTr', ...
        'omzICpostTr', 'zetaICpostTr','-v7.3');
end

end

