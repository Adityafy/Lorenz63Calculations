function u = stdynamics(p,u,solverfunc,linpartfunc,nonlinpartfunc)
% u = stdynamics(p,u,solverfunc,linpartfunc,nonlinpartfunc)
% calculates the space time dynamics (attractor)
linpart = linpartfunc(p);
nmax = p.nmax;
for n = 1:nmax
    u(:,n+1) = solverfunc(p,u(:,n),n,linpart,nonlinpartfunc);
end
end