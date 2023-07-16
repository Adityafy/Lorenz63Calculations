function [v,R,laminst,lamgs] = gs(para,h,nnorm,dynamics,evolFunc,jacobianFunc)
% [v,R,laminst,lamgs] = gs(paramVec,dynamicsVec,coupMatFunc,jacobianFunc)
% feed one time-step dynamics to the jacobian
[M,total_time] = size(dynamics);
v(:,:,1) = eye(M); % initial perturbation vectors
R(:,:,1) = zeros(M,M);
laminst = zeros(M,1);
lamgs = zeros(M,1);
% coupling_matrix = coupMatFunc(paramVec);
fprintf("Calculating GS----------\n")
for t = 1:total_time
    % J = jacobianFunc(para,dynamics(:,t));
    pertVecs = evolFunc(para,h,dynamics(:,t),v(:,:,t),jacobianFunc);
    [v(:,:,t+1), R(:,:,t)] = qr(pertVecs);
    if rem(t, nnorm) == 0
        for k = 1:M
            laminst(k,t/nnorm) = log(abs(R(k,k,t)));
        end
    end
    if rem(t, total_time/100) == 0
        fprintf("This is time step %g and ", t);
        toc
    end
end
toc
for k = 1:M
    lamgs(k,:) = (nnorm/total_time)*sum(laminst(k,:));
end

end