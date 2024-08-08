function omztilde = rk2gsh2(omz,psi,p,dynfunc, h)
% omztilde = rk2gsh2(omz,psi,p,dynfunc,run_type)
% run_type determines the type of run, whether it is a run for transients
% (run_type = 0) or if it's for the main dynamics (run_type = 1).
    % if run_type == 0
    %     h = p.dt_tr;
    % elseif run_type == 1
    %     h = p.dt;
    % end
    % h = p.dt;
    k1 = dynfunc(omz,psi,p);
    %u1 = u+k1*dt;
    k2 = dynfunc(omz+k1*h,psi,p);
    omztilde = omz + h*((k1+k2)/2);
end