function cmat = cmatrix(R,checkCond,backSub)
%
% cmat = cmatrix(p,R,checkCond,backSub)
% Calculates C matrix for the calculation of CLVs. 
% p = parameters
% R = R matrix
% checkCond = 0 or 1 : switch to calculate the condition number of
%                       the R matrix
% backSub = 0 or 1 : switch to calculate c matrix by doing a back
%                       subsitution instead of matrix inversion

M = size(R,2);
tu = size(R,3);
cnorm = 1; % normalization time for backward evolution
rng(2);
cmat = zeros(M,M,tu+1);
for i = 1:M
    cmat(1:i,i,tu+1) = rand(i,1);
end
fprintf("Calculating c matrix----------\n")
if backSub == 0
    for j = tu:-1:1
        for i = 1:M
            if checkCond == 1
                condNum = cond(R(1:i,1:i,j));
                % if det(R(1:i,1:i,j)) == 0
                %     error('det(R(1:i,1:i,j)) is zero');
                % end
                if condNum > 1e4
                    fprintf(['Time step: %d, Matrix size: %d, ' ...
                        'High condition number: %e \n'], j, i, condNum);
                end
            end
            cmat(1:i,i,j) = R(1:i,1:i,j)\cmat(1:i,i,j+1);
        end
        if j~=tu && rem(j,cnorm) == 0
            for i = 1:M
                cmat(1:i,i,j) = cmat(1:i,i,j)./norm(cmat(1:i,i,j));
            end
        end
        if rem(j, tu/10) == 0
            fprintf("This is time step %g and ", j);
            toc
        end
    end
elseif backSub == 1
    for t = tu:-1:1
        for n = 1:M
            cmat(n,n,t) = cmat(n,n,t+1)/R(n,n,t);
            for j = 1:n-1
                sum = 0;
                for i = n-(j-1):n
                    sum = sum + R(n-j,i,t) * cmat(i,n,t);
                end
                cmat(n-j,n,t) = (cmat(n-j,n,t+1)-sum) / R(n-j,n-j,t);
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