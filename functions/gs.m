function [v,R,laminst,lamgs] = gs(p,dynamics,solverfunc,jacobianFunc)
% [v,R,laminst,lamgs] = gs(p,dynamics,solverfunc,jacobianFunc)
% Performs the Gram-Schmidt reorthonormalization.
% Calculates the evolution of perturbation vectors v, the R matrix 
% from QR, and the instantaneous and time-averaged Lyapunov exponents.

nnorm1 = p.nnorm1;
tnorm = p.tnorm;

[M,nmax] = size(dynamics); % M is system dimension
v(:,:,1) = eye(M); % initial perturbation vectors
R(:,:,1) = zeros(M,M);
laminst = zeros(M,1);
lamgs = zeros(M,1);

nonlinpart = @(v,n) [0; 0; 0];

fprintf("Calculating GS----------\n")
for n = 1:nmax
    linpart = jacobianFunc(p,dynamics(:,n));
    pertvecs = solverfunc(p,v(:,:,n),n,linpart,nonlinpart);
    [v(:,:,n+1), R(:,:,n)] = qr(pertvecs);
    if rem(n, nnorm1) == 0
        for k = 1:M
            laminst(k,n/nnorm1) = (1/tnorm)*log(abs(R(k,k,n)));
        end
    end
end

for k = 1:M
    lamgs(k,:) = (1/nmax)*sum(laminst(k,:));
end

end