function cmat = cmatrix(p,x,R,checkCond,backSub)
%
% cmat = cmatrix(x,R)
% Calculates C matrix for the calculation of CLVs. Requires the lattice
% dynamics x and R matrix as input. x is a matrix with number of rows as 
% system size and number of columns as time-steps
[M,total_time] = size(x); % M = Nx*Ny
cnorm = p.cnorm; % normalization time for backward evolution

rng(2);
cmat = zeros(M,M,total_time+1);
for i = 1:M
    cmat(1:i,i,total_time+1) = rand(i,1);
end
fprintf("Calculating c matrix----------\n")
if backSub == 0
    for t = total_time:-1:1
        for i = 1:M
            if checkCond == 1
                condNum = cond(R(1:i,1:i,t));
                % if det(R(1:i,1:i,j)) == 0
                %     error('det(R(1:i,1:i,j)) is zero');
                % end
                if condNum > 1e4
                    fprintf('Time step: %d, Matrix size: %d, High condition number: %e \n', t, i, condNum);
                end
            end
            cmat(1:i,i,t) = R(1:i,1:i,t)\cmat(1:i,i,t+1);
        end
        if t~=total_time && rem(t,cnorm) == 0
            for j = 1:M
                cmat(:,j,t) = cmat(:,j,t)./norm(cmat(:,j,t));
            end
        end
        if rem(t, total_time/10) == 0
            fprintf("This is time step %g and ", t);
            toc
        end
    end
elseif backSub == 1
    for t = total_time:-1:1
        for n = 1:M
            cmat(n,n,t) = cmat(n,n,t+1)/R(n,n,t);
            for t = 1:n-1
                sum = 0;
                for i = n-(t-1):n
                    sum = sum + R(n-t,i,t) * cmat(i,n,t);
                end
                cmat(n-t,n,t) = (cmat(n-t,n,t+1)-sum) / R(n-t,n-t,t);
            end
        end
        if t~=total_time && rem(t,cnorm) == 0
            for i = 1:M
                cmat(1:i,i,t) = cmat(1:i,i,t)./norm(cmat(1:i,i,t));
            end
        end
        if rem(t, total_time/10) == 0
            fprintf("This is time step %g and ", t);
            toc
        end
    end
end
toc
end