function p = paramgsh(eps,sig,csq,gm,lam_0,dx,dy,dt_tr,dt_fpv,gamma, ...
    trtimeu,totimeu,tN,seed)
% p = paramgsh(eps,sig,csq,gm,dx,dy,dt,rolls,trtimeu,totimeu,seed)
% lam_0 is 2*pi
% dx or dy should be taken as 2*pi/k, k is an integer


c = sqrt(csq);
totmax = totimeu/dt_fpv; % total time steps for dynamics after transient

% ----domain size----
% Nx = round(rolls*lam_0/dx); % spatial nodes in x
% Ny = round(rolls*lam_0/dy); % spatial nodes in y
Nx = round(0.5*gamma*lam_0/dx);
Ny = round(0.5*gamma*lam_0/dy);

N = Nx*Ny;

gmPos = 2; % position of gm; 1 for eqn 1, 2 for eqn 2 (Karimi)
gsiters = 5; % number of iterations for gauss seidel

% ----random seed----
rng(seed,"twister");

% ------- running filename ---------
run_name = join([num2str(Nx) 'x' num2str(Ny) 'eps' dashed(eps) ...
            'sig' dashed(sig) 'csq' dashed(csq) ...
            'gm' dashed(gm) 'tN' dashed(tN) 'trt' dashed(trtimeu) ...
            'tot' dashed(totimeu) 's' dashed(seed)]);

% ----video----
altframes = 2;

% ------ vectors for periodic BCs -----------
I = 1:Nx;
J = 1:Ny;
Ip1 = circshift(I,-1);
Im1 = circshift(I,1);
Ip2 = circshift(I,-2);
Im2 = circshift(I,2);
Jp1 = circshift(J,-1);
Jm1 = circshift(J,1);
Jp2 = circshift(J,-2);
Jm2 = circshift(J,2);

%-------- intervals ----------
% if totmax < intervals
%     error('Number of intervals is less than total time steps!');
%     fprintf('Total time steps = %g', p.tmax);
%     fprintf('Number of intervals = %g', p.intervals);
% end

TrIntervals = floor(trtimeu);
ToIntervals = floor(totimeu);
if totimeu < 1
    ToIntervals = 10;
end


% ----  (p)arameter (struct)ure ----
p = struct('eps', eps, 'sig', sig, 'c', c, 'gm', gm, 'dt_tr', dt_tr, ...
    'dt_fpv', dt_fpv, 'dx', dx, 'dy', dy, 'tmax', totmax, 'Nx', Nx, ...
    'Ny', Ny, 'N',N, 'seed', seed, 'I', I, 'J', J, 'Ip1', Ip1, ...
    'Im1', Im1,'Ip2', Ip2, 'Im2', Im2,'Jp1', Jp1,'Jm1', Jm1,'Jp2', Jp2, ...
    'Jm2', Jm2, 'gmPos', gmPos, 'run_name', run_name, 'lam_0', lam_0, ...
    'trtimeu', trtimeu, 'totimeu', totimeu, 'tN', tN, 'gamma', gamma, ...
    'TrIntervals', TrIntervals,'ToIntervals',ToIntervals);

    function outstr = dashed(parameter)
        parastring = num2str(parameter);
        if contains(parastring,'.')
            strarr = split(parastring,'.');
            outstr = join(strarr,'-');
            outstr = char(outstr);
        else
            outstr = char(parastring);
        end
    end

end

