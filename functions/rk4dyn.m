function dynamics = rk4dyn(p,dynamics,dynrhsfunc)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
dt = p.dt;
nmax = p.nmax;
for j = 1:nmax
    % k values
    k1 = dynrhsfunc(p, dynamics(:,j));
    k2 = dynrhsfunc(p, dynamics(:,j) + (0.5*dt)*k1);
    k3 = dynrhsfunc(p, dynamics(:,j) + (0.5*dt)*k2);
    k4 = dynrhsfunc(p, dynamics(:,j) + dt*k3);
    % dynamics
    dynamics(:,j+1) = dynamics(:,j) + (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);
end
end

