function [psimat, omzmat, dH1, dH1mag,res] = advGSHdH1stepFDSI( ...
    psimat, omzmat, zetamat, dpsi1mat, domz1mat, dzeta1mat, ...
    matdivpsi, matdivomz, n, dH1mag, res_calc, res, p)

%ADVGSHDH1STEPFDSI Advances the state of GSH and tangent space using RK2 and Crank-Nicolson.
%
%   This function performs one step of:
%   - Updating the zeta field iteratively from vorticity
%   - Computing nonlinear evolution using RK2
%   - Applying implicit Crank-Nicolson step for linear part
%   - Advancing tangent space equations with optional residual calculation
%
%   Inputs:
%     psimat, omzmat      : Current state fields
%     zetamat             : Initial zeta guess
%     dpsi1mat, domz1mat  : Perturbation fields
%     dzeta1mat           : Initial perturbation zeta guess
%     matdivpsi, matdivomz: Crank-Nicolson inversion matrices
%     dt                  : Time step
%     Nx, Ny              : Grid dimensions
%     n                   : Current time index (for fpvmag)
%     fpvmag              : Array storing magnitudes of perturbation vectors
%     res_calc            : Flag (0 = no residual, 1 = compute residual)
%     res                 : Residual array (used if res_calc == 1)
%     p                   : Parameter struct
%
%   Outputs:
%     psimat, omzmat      : Updated state fields
%     dpsi1mat, domz1mat  : Updated perturbation fields
%     fpvmag              : Updated perturbation magnitude array
%     res                 : Updated residual array

dt = p.sim.dt;
Nx = p.rmesh.Nx;
Ny = p.rmesh.Ny;

%--- iterative inversion of zeta from vorticity
zetamat = iterativeZetaOmz(omzmat, zetamat, p);
% zetamat = zetaGSHspectral(p,omzmat);
% [ulat, vlat] = uvzeta(zetamat, p);

%--- nonlinear evolution via RK2
psitilde = rk2gsh1(psimat, zetamat, p, @nlpgsh1, dt);
psitilde = reshape(psitilde',[],1);

omztilde = rk2gsh2(omzmat, psimat, p, @nlpgsh2, dt);
omztilde = reshape(omztilde',[],1);

%--- Crank-Nicolson implicit step for linear part
psivec = matdivpsi * psitilde;
omzvec = matdivomz * omztilde;

psimat = reshape(psivec, Nx, Ny)';
omzmat = reshape(omzvec, Nx, Ny)';

%--- tangent space dynamics
dzeta1mat = iterativeZetaOmz(domz1mat, dzeta1mat, p);

if res_calc == 0
    [dpsi1mat, ~] = rk2tsgsh1(psimat, zetamat, dpsi1mat, ...
        dzeta1mat, p, @ts1exp, @sfdd);
elseif res_calc == 1
    [dpsi1mat, res(n)] = rk2tsgsh1(psimat, zetamat, dpsi1mat, ...
        dzeta1mat, p, @ts1exp, @sfdd);
end

domz1mat = rk2tsgsh2(omzmat, psimat, domz1mat, ...
    dpsi1mat, p, @ts2exp, @sfdd);

%--- norm of full perturbation vector
dH1 = [reshape(dpsi1mat',[],1); reshape(domz1mat',[],1)];
dH1mag(n+1) = norm(dH1);
end
