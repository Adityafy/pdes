function [clv,cle,clvmag] = computeCLVs(p,gsVec,R,checkCond,backSub)

tu = size(R,3);
c_conv_time = round(0.4*tu); % time for random c coeffs to converge
clv_time = tu - c_conv_time; % Define clvtime based on c_conv_time

nv = p.ts.nv;
N = p.rmesh.N;
% Nx = p.rmesh.Nx;
% Ny = p.rmesh.Ny;
clv = zeros(2*N,nv,clv_time); % w (clv) initial

M = size(R,2);

cnorm = 1; % normalization time for backward evolution
rng(2);
cmat = zeros(M,M,tu+1);
clvmag = zeros(nv,1);

for i = 1:M
    cmat(1:i,i,tu+1) = rand(i,1);
end
fprintf("Going back in time, calculating c matrix... \n")
if backSub == 0
    % backward time evolution
    for n = tu:-1:1
        for i = 1:M
            if checkCond == 1
                condNum = cond(R(1:i,1:i,n));
                % if det(R(1:i,1:i,j)) == 0
                %     error('det(R(1:i,1:i,j)) is zero');
                % end
                if condNum > 1e4
                    fprintf(['Time step: %d, Matrix size: %d, ' ...
                        'High condition number: %e \n'], n, i, condNum);
                end
            end
            cmat(1:i,i,n) = R(1:i,1:i,n)\cmat(1:i,i,n+1);
        end
        if n~=tu && rem(n,cnorm) == 0
            for i = 1:M
                cmat(1:i,i,n) = cmat(1:i,i,n)./norm(cmat(1:i,i,n));
            end
        end
        % check if n is less than or equal to clv_time then start
        % calculating CLVs
        if n <= clv_time
            if n == clv_time
                fprintf("Backward time evolution continues,\nNow calculating cmatrix and CLVs");
            end
            % for every nv calculate the clv as a linear comb of gsVec with
            % cmat coefficients
            for j  = 1:nv
                for k = 1:j
                    clv(:,j,n) = clv(:,j,n) + cmat(k,j,n) .* gsVec(:,k,n);
                end
            end
            % calculate cle and renormalize CLVs periodically
            if rem(n,p.ts.nnorm) == 0
                for k = 1:nv
                    clv_mag_n(k) = norm(clv(:,k,n));
                    cle(k,n) = (1/p.ts.tN) * log(abs(clv_mag_n(k)));
                    clv(:,k,n) = clv(:,k,n)./clv_mag_n(k);
                end
                clvmag(:,n) = clv_mag_n(k);
            end
        end
        % time display
        if rem(n, tu/10) == 0
            fprintf("This is time step %g and ", n);
            toc
        end
    end
elseif backSub == 1
    for t = tu:-1:1
        for dim = 1:M
            cmat(dim,dim,t) = cmat(dim,dim,t+1)/R(dim,dim,t);
            for dim = 1:dim-1
                sum = 0;
                for i = dim-(dim-1):dim
                    sum = sum + R(dim-dim,i,t) * cmat(i,dim,t);
                end
                cmat(dim-dim,dim,t) = (cmat(dim-dim,dim,t+1)-sum) / R(dim-dim,dim-dim,t);
            end
        end
        if t~=tu && rem(t,cnorm) == 0
            for i = 1:M
                cmat(1:i,i,t) = cmat(1:i,i,t)./norm(cmat(1:i,i,t));
            end
        end
        if rem(t, tu/10) == 0
            fprintf("This is time step %g and ", t);
            toc
        end
    end
end
toc
end