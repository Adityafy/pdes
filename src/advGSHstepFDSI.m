function [psimat, omzmat, zetamat, ulat, vlat] = advGSHstepFDSI(psimat, omzmat, ...
                                         zetamat, matdivpsi, matdivomz, p)
% advGSHstepFDSI Advance one time step of the GSH equation using Finite
% Difference and Semi-Implicit (Crank Nicolson)
%
%   [psimat, omzmat, ulat, vlat] = advGSHstepFDSI(psimat, omzmat, 
%       matdivpsi, matdivomz, p, dt)
%
%   Inputs:
%       psimat     - current psi field (Nx x Ny)
%       omzmat     - current omega_z field (Nx x Ny)
%       matdivpsi  - implicit Crank-Nicolson matrix for psi
%       matdivomz  - implicit Crank-Nicolson matrix for omega_z
%       p          - parameter struct with fields like Nx, Ny, etc.
%       dt         - time step size
%
%   Outputs:
%       psimat     - updated psi field
%       omzmat     - updated omega_z field
%       ulat, vlat - velocity field (from streamfunction)
%
%   Author: Aditya

    % -- Zeta update (iterative) --
    zetamat = iterativeZetaOmz(omzmat, zetamat, p);
    % zetamat = zetaGSHspectral(p,omzmat);
    [ulat, vlat] = uvzeta(zetamat, p);

    % -- Nonlinear explicit RK2 update --
    psitilde = rk2gsh1(psimat, zetamat, p, @nlpgsh1, p.sim.dt);
    psitilde = reshape(psitilde', [], 1);  % latToVec replacement

    omztilde = rk2gsh2(omzmat, psimat, p, @nlpgsh2, p.sim.dt);
    omztilde = reshape(omztilde', [], 1);

    % -- Linear implicit Crank-Nicolson --
    psivec = matdivpsi * psitilde;
    omzvec = matdivomz * omztilde;

    % -- Convert back to matrix form --
    psimat = reshape(psivec, p.rmesh.Ny, p.rmesh.Nx)';
    omzmat = reshape(omzvec, p.rmesh.Ny, p.rmesh.Nx)';
end
