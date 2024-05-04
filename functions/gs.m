function [v,R,laminst,lamgs] = gs(p,dynamics,evolFunc,jacobianFunc)
% [v,R,laminst,lamgs] = gs(p,dynamics,evolFunc,jacobianFunc)
% feed one time-step dynamics to the jacobian

nnorm1 = p.nnorm1;
tnorm = p.tnorm;

[M,nmax] = size(dynamics);
v(:,:,1) = eye(M); % initial perturbation vectors
R(:,:,1) = zeros(M,M);
laminst = zeros(M,1);
lamgs = zeros(M,1);

fprintf("Calculating GS----------\n")
for t = 1:nmax
    pertvecs = evolFunc(p,dynamics(:,t),v(:,:,t),jacobianFunc);
    [v(:,:,t+1), R(:,:,t)] = qr(pertvecs);
    if rem(t, nnorm1) == 0
        for k = 1:M
            laminst(k,t/nnorm1) = (1/tnorm)*log(abs(R(k,k,t)));
        end
    end
end

for k = 1:M
    lamgs(k,:) = (1/nmax)*sum(laminst(k,:));
end

end