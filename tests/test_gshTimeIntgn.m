%% make sure to run ../main/gsh.m with the same parameters
addpath('../main/');
mfa = '/media'; % media folder address
dfa = '/data'; % saved data folder address

tic
run('../main/gsh.m');
psimat_oldend = psimat(:,:,end);

%%

% adjusted parameters for new run to make it same as gsh

trtimeu = round(time_units/2); % transient time units
totimeu = trtimeu; % total time units


p = paramgsh(eps,sig,csq,gm,dx,dy,dt,rolls,trtimeu,totimeu,seed);

%% 

[psimat, omzmat, zetamat, ulat, vlat] = gshTimeIntgn(p,@gshtric);
toc

psimat_newend = psimat(:,:,end);


%%
testpsi = psimat_oldend == psimat_newend;
if sum(sum(testpsi)) == Nx*Ny
    sum(sum(testpsi))
    fprintf('\nPASSED!!\n')
end