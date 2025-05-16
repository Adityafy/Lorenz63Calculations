clear all;
addpath('../functions/');

%% parameters
sigma = 10; r = 28; b = 8/3;
%lorPara = [sigma, r, b]; % lorenz parameters

time_u = 500;

dt = 0.01;
nmax = time_u/dt;
timesteps_vec = 0:dt:time_u;
h = dt;

ic = [-7.1778; -12.9972; 12.5330]; % initial conditions
u = zeros(max(size(ic)), nmax); % state vector
u(:,1) = ic;

[m,nmax] = size(u); % m is the system dimension
M = m*nmax; 


% Gram-Schmidt parameters
nnorm1 = 1;
tnorm = nnorm1*dt;

% c calculation parameters
nnorm2 = 10;

p = struct('sigma', sigma, 'r', r, 'b', b, 'nmax', nmax, 'dt', dt, 'm', m, ...
            'nnorm1', nnorm1, 'nnorm2', nnorm2, 'tnorm', tnorm, 'time_u', time_u, ...
            'timesteps_vec', timesteps_vec);

%% exponential time integration (Cox and Matthews)
u = etdrk4(p,u,@linLorenz63,@nonlinLorenz63);
% u = etd1(p,u,@linLorenz63,@nonlinLorenz63);

%% dynamics plot
fig_attractor(u);

% %% ETD cross1994
% function u = etdcross94(p,u,Lfunc,Nfunc)
% dt = p.dt;
% m = p.m;
% nmax = p.nmax;
% t = p.timesteps_vec;
% L = Lfunc(p);
% N = Nfunc;
% 
% for n = 1:nmax-1
%     un = u(:,n);
%     u(:,n+1) = expm(L*dt)*un + N(un)*(expm(L*dt)-eye(m))
% end
% end

%% ETD1 function
function u = etd1(p,u,Lfunc,Nfunc)

dt = p.dt;
m = p.m;
nmax = p.nmax;
t = p.timesteps_vec;
L = Lfunc(p);
N = Nfunc;

invL = inv(L);

for n = 1:nmax-1
    un =  u(:,n);
    u(:,n+1) = expm(L*dt)*un + invL * ((expm(L*dt)-eye(m)) * N(un,t(n)));
end

end

%% ETD RK4 function
function u = etdrk4(p,u,Lfunc,Nfunc)
% u = etdrk4(p,u,Lfunc,Nfunc)

h = p.dt;
m = p.m;
nmax = p.nmax;
t = p.timesteps_vec;

% precomputing some
L = Lfunc(p);
N = Nfunc;
invL = inv(L);

E1 = expm(L*h);
E2 = expm(L*h/2);
omega = invL * (E2-eye(m));

eta = h^(-2) * L^(-3);
alpha = eta * (- 4*eye(m) - L*h + E1 * (4*eye(m) - 3*L*h + (L*h)^2));
beta = eta * (2*eye(m) + L*h + E1 * (-2*eye(m) + L*h));
gamma = eta * (-4*eye(m) - 3*L*h - (L*h)^2 + E1 * (4*eye(m) - L*h));

% time integration
for n = 1:nmax
    an = E2 * u(:,n) + omega * N(u(:,n),t(n));
    bn = E2 * u(:,n) + omega * N(an,t(n)+h/2);
    cn = E2 * an + omega * ( 2 * N(bn,t(n)+h/2) - N(u(:,n),t(n)) );

    u(:,n+1) =  E1 * u(:,n) + alpha * N(u(:,n), t(n)) ...
        + beta * 2 * ( N(an, t(n)+h/2) + N(bn, t(n)+h/2) ) ...
        + gamma * N(cn, t(n)+h) ;

end

end