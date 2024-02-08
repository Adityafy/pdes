
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
dt = 0.2;
rolls = 5;
trtimeu = 50; % transient time units
totimeu = 50; % total time units
seed = 1;

p = paramgsh(eps,sig,csq,gm,dx,dy,dt,rolls,trtimeu,totimeu,seed);

%% Dynamics from initial state
tic
[psimat, omzmat, zetamat, ulat, vlat] = gshTimeIntgn(p,@gshtric);
toc

%%
run('postprocessgsh.m');