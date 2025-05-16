function x = rk4dyn(p,x,dynrhsfunc,linpart,nonlinpart)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
dt = p.dt;
nmax = p.nmax;
for j = 1:nmax
    % k values
    k1 = dynrhsfunc(p, x(:,j) , linpart, nonlinpart);
    k2 = dynrhsfunc(p, x(:,j) + (0.5*dt)*k1 , linpart, nonlinpart);
    k3 = dynrhsfunc(p, x(:,j) + (0.5*dt)*k2 , linpart, nonlinpart);
    k4 = dynrhsfunc(p, x(:,j) + dt*k3 , linpart, nonlinpart);
    % dynamics
    x(:,j+1) = x(:,j) + (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);
end
end

