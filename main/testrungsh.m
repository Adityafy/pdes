
clear all;

addpath('../src/');

%% control parameters
epsilon = 0.7;
sig = 1;
csq = 0.5;
gm = 50;

% parameters for spectral (ps)
cp = struct('epsilon',epsilon,'sig',sig,'csq',csq,'gm',gm);

%% spatial variables
Lx = 64;
sv = spatialcalcsgshspec(Lx,cp);

%% transient parameters (trp)
dt = 0.1;
trtu = 100;
seed = 1;
trp = struct('dt',dt,'trtu',trtu,'seed',seed);

run_name = join(['specTr' num2str(sv.Nx) 'eps' pointToDashed(cp.epsilon) ...
            'sig' pointToDashed(cp.sig) 'csq' pointToDashed(cp.csq) ...
            'gm' pointToDashed(cp.gm) 'trt' pointToDashed(trp.trtu) ...
            's' pointToDashed(trp.seed)]);

%% GSH spectral transients
[psitr,omztr,zetatr,utr,vtr] = gshSpecETDTrFull(cp,sv,trp);