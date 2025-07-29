function etd = etdPreCalcs(L1,L2,dt)
% [expL1dtmat,expL2dtmat,expL1dtvec,expL2dtvec] = etdprecalcs(L1,L2,dt)
% Calculates the static matrices used in ETD1
expL1dtmat = exp(L1*dt);
expL2dtmat = exp(L2*dt);
expL1dtvec = exp(reshape(L1',[],1)*dt);
expL2dtvec = exp(reshape(L2',[],1)*dt);
etd = struct('expL1dtmat',expL1dtmat,'expL2dtmat',expL2dtmat, ...
                'expL1dtvec',expL1dtvec,'expL2dtvec',expL2dtvec);
end