function psitilde = rk2gsh1(psi,zeta,p,dynfunc,h)
% psitilde = rk2gsh1(psi,zeta,p,dynfunc,run_type)
% run_type determines the type of run, whether it is a run for transients
% (run_type = 0) or if it's for the main dynamics (run_type = 1).
    % if run_type == 0
    %     h = p.dt_tr;
    % elseif run_type == 1
    %     h = p.dt;
    % end

    % h = p.dt;
    k1 = dynfunc(psi,zeta,p);
    %u1 = u+k1*dt;
    k2 = dynfunc(psi+k1*h,zeta,p);
    psitilde = psi + h*((k1+k2)/2);
end