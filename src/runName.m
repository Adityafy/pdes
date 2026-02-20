function dyn_run_name = runName(p)

if p.sim.runtype == 0
    runstr = 'fd';
elseif p.sim.runtype == 1
    runstr = 'ps';
else
    error('Proper runtype not specified (simp.runtype = 0 for FD, 1 for PS).');
end

dyn_run_name = join([runstr 'G' num2str(p.sim.Gamma) 'eps' pointToDashed(p.con.epsilon) ...
            'sig' pointToDashed(p.con.sig) 'csq' pointToDashed(p.con.csq) ...
            'gm' pointToDashed(p.con.gm) 'time' pointToDashed(p.sim.tstart) '-' ...
            pointToDashed(p.sim.tend) 's' pointToDashed(p.sim.seed)]);

end