function u = rk4(p,u,n,linpart,nonlinpartfunc)
% u = rk4(p,u,n,linpart,nonlinpartfunc)
% calculates one time step for rk4
% requires a linear part (operator) and 
% nonlinear part function

dynrhsfunc = @(p,x,linpart,nonlinpartfunc) linpart*x+nonlinpartfunc(x,n);
dt = p.dt;

% k values
k1 = dynrhsfunc(p, u , linpart, nonlinpartfunc);
k2 = dynrhsfunc(p, u + (0.5*dt)*k1 , linpart, nonlinpartfunc);
k3 = dynrhsfunc(p, u + (0.5*dt)*k2 , linpart, nonlinpartfunc);
k4 = dynrhsfunc(p, u + dt*k3 , linpart, nonlinpartfunc);

% dynamics
u = u + (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);

end