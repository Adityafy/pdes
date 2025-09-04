function [psiICpostTr, omzICpostTr, zetaICpostTr] = gshTimeIntgICpostTr(p, ICfunc)
% gshTimeIntgICpostTr
% ------------------------------------------------------------
% Time-integrates the Generalized Swift-Hohenberg (GSH) system 
% for a fixed number of steps (transient removal), and outputs 
% the evolved fields as new initial conditions.
%
% [psiICpostTr, omzICpostTr, zetaICpostTr, ...
%    uICpostTr, vICpostTr] = gshTimeIntgICpostTr(p, ICfunc)
%
% Inputs:
%   p       - struct containing simulation parameters and settings
%   ICfunc  - function handle to generate initial conditions
%
% Outputs:
%   psiICpostTr   - psi field after transients
%   omzICpostTr   - omega_z field after transients
%   zetaICpostTr  - zeta field after transients
%   uICpostTr     - u-velocity field after transients (unused here, returned as empty)
%   vICpostTr     - v-velocity field after transients (unused here, returned as empty)
%
% Author: Aditya

% Total number of time steps to integrate (transient duration)
nmax  = p.sim.nmax;

% If using the finite-difference semi-implicit scheme, precompute 
% implicit matrix inverses to accelerate time-stepping.
if p.sim.runtype == 0
    [matdivpsi, matdivomz] = impMatGSH(p, p.sim.dt);
end

% Initialize simulation state using provided IC function
% Only psi, omz, and zeta fields are needed for transient integration
[psimat, ~, omzmat, ~, zetamat, ~, ~, ~] = ICfunc(p, p.sim.dynICtype);

% for etd
etd = etdPreCalcs(p.L1,p.L2,p.sim.dt);

% =================== Begin Time Integration ===================
for n = 1:nmax
    if p.sim.runtype == 0
        % -- Finite Difference Semi-Implicit Scheme --
        % Advance fields using finite-difference time stepping
        [psimat, omzmat, zetamat] = advGSHstepFDSI( ...
            psimat, omzmat, zetamat, matdivpsi, matdivomz, p);

    elseif p.sim.runtype == 1
        % -- Pseudospectral Exponential Time Differencing Scheme --
        
        % Compute zeta from omega_z in Fourier space (Poisson-like solve)
        zetamat = zetaGSHspectral(p, omzmat);

        % Advance fields using ETD-based pseudospectral method
        [psimat, omzmat, zetamat] = advGSHstepPSETD( ...
            p, psimat, omzmat, zetamat, ...
            etd.expL1dtmat, etd.expL2dtmat, @N1hat, @N2hat);
    end

    % ========== Optional Progress Report ==========
    if rem(n, round(nmax / p.sim.progReportFactor)) == 0
        fprintf('This is time step: %g / %g, ', n, nmax);
        toc;
    end

    % ========== NaN Check for Numerical Blow-Up ==========
    if isnan(abs(sum(sum(psimat))))
        error('Blow up occurred in psi..');
    end

    % ========== Optional Live Visualization ==========
    if p.sim.makeLiveFig == 1
        imagesc(psimat); clim([-1 1]);          % display psi field
        % imagesc(zetamat); clim([-10 10]);
        colorbar;               % add color scale
        axis square;            % force square aspect ratio
        %          % fix color scale for visibility
        colormap jet;           % set colormap
        title(join(['t = ', num2str(p.sim.dt*n)]));
        drawnow;                % update figure window
    end
end
% =================== End of Time Integration ===================

% Assign final state after transients as "new" initial conditions
psiICpostTr  = psimat;
omzICpostTr  = omzmat;
zetaICpostTr = zetamat;
% uICpostTr    = umat;   % Not evolved or used in this function
% vICpostTr    = vmat;   % Not evolved or used in this function
end
