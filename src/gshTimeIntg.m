function [psi, omz, zeta, u, v] = gshTimeIntg(p, ICfunc)
%GSHTIMEINTG Integrates the Generalized Swift-Hohenberg (GSH) system in time.
%
%   [psi, omz, zeta, u, v] = gshTimeIntg(p, ICfunc)
%
%   Inputs:
%     p       : Struct containing simulation parameters and numerical options
%     ICfunc  : Function handle to generate initial conditions
%
%   Outputs:
%     psi     : Saved snapshots of psi at output intervals
%     omz     : Saved snapshots of vorticity at output intervals
%     zeta    : Saved snapshots of streamfunction at output intervals
%     u, v    : Saved snapshots of velocity components at output intervals
%
%   Author: Aditya

    nmax  = p.sim.nmax;     % Total number of time steps
    interv = p.sim.interv;  % Number of saved output intervals
    etd = etdPreCalcs(p.L1,p.L2,p.sim.dt); % for etd

    % Error if time steps are too large
    if nmax >= 100000
        error("Time steps are too large, this transient run is not " + ...
            "recommended. Use the transient run that doesn't save transients.");
    end

    % Precompute implicit matrix inverses if using finite-difference method
    if p.sim.runtype == 0
        [matdivpsi, matdivomz] = impMatGSH(p, p.sim.dt);
    end

    % Initialize state and diagnostics using provided initial condition function
    [psimat, psi, omzmat, omz, zetamat, zeta, u, v] = ICfunc(p, p.sim.dynICtype);

    for n = 1:nmax
        if p.sim.runtype == 0
            % Advance using finite-difference semi-implicit scheme
            [psimat, omzmat, zetamat, umat, vmat] = advGSHstepFDSI( ...
                psimat, omzmat, zetamat, matdivpsi, matdivomz, p);
        elseif p.sim.runtype == 1
            % Compute zeta using pseudospectral method
            zetamat = zetaGSHspectral(p, omzmat);

            % Advance using pseudospectral exponential time differencing (ETD)
            [psimat, omzmat, zetamat, umat, vmat] = advGSHstepPSETD( ...
                p, psimat, omzmat, zetamat, ...
                etd.expL1dtmat, etd.expL2dtmat, @N1hat, @N2hat);
        end

        % ============== Save simulation snapshots at specified intervals ==============
        if rem(n, nmax / interv) == 0
            indx = round(n * interv / nmax);
            psi(:, :, indx)  = psimat;
            omz(:, :, indx)  = omzmat;
            zeta(:, :, indx) = zetamat;
            u(:, :, indx)    = umat;
            v(:, :, indx)    = vmat;
        end

        % ============== Progress update every 10% of simulation ==============
        if rem(n, round(nmax / p.sim.progReportFactor)) == 0
            fprintf('This is time step: %g / %g, ', n, nmax);
            toc;
        end

        % ============== Check for numerical blow-up ==============
        if isnan(abs(sum(sum(psimat))))
            error('Blow up occurred in psi..');
        end

        % ============== Optional live figure ==============
        if p.sim.makeLiveFig == 1
            imagesc(psimat);
            colorbar;
            axis square;
            clim([-1 1]);
            colormap jet;
            drawnow;
        end
    end
end
