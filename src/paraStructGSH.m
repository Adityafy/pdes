function p = paraStructGSH(con,sim,spatialGridRealFunc,...
                            spatialGridSpecFunc,LinOpFunc)
% p = paraStructGSH(con,sim,spatialGridRealFunc,...
%                           spatialGridSpecFunc,LinOpFunc)
p.con = con;
p.sim = sim;
p.rmesh = spatialGridRealFunc(p.sim.lam_0, p.sim.wlmult);
p.smesh = spatialGridSpecFunc(p.rmesh.Lx, p.rmesh.Nx);
p.idx = makeIndexStructPBC(p.rmesh.Nx,p.rmesh.Ny); % for PBC indexing in FDSI
[p.L1,p.L2] = LinOpFunc(p,p.smesh.Kx, p.smesh.Ky);

if sim.runtype == 0
    runstr = 'fd';
elseif sim.runtype == 1
    runstr = 'ps';
else
    error('Proper runtype not specified (simp.runtype = 0 for FD, 1 for PS).');
end

dyn_run_name = join([runstr num2str(ceil(p.rmesh.Lx)) 'N' num2str(p.rmesh.Nx) 'eps' pointToDashed(p.con.epsilon) ...
            'sig' pointToDashed(p.con.sig) 'csq' pointToDashed(p.con.csq) ...
            'gm' pointToDashed(p.con.gm) 'tu' pointToDashed(p.sim.tu) ...
            's' pointToDashed(p.sim.seed)]);

p.dyn_run_name = dyn_run_name;

p.Kdiff = (p.smesh.Kx).^2 + (p.smesh.Ky).^2;
firstrowKdiff = latToVec(p.Kdiff)';
p.KdiffMat = toeplitz(firstrowKdiff);
p.qx = latToVec(p.smesh.Kx);
p.qy = latToVec(p.smesh.Ky);

p.etd = etdPreCalcs(p.L1,p.L2,p.sim.dt);

end