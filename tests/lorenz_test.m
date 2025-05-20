% Lorenz CLV calculations runscript
close all;
clear all;
% addpath('functions/');

%% parameters
sigma = 10; r = 28; b = 8/3;
%lorPara = [sigma, r, b]; % lorenz parameters

time_u = 200;

dt = 0.001;
nmax = time_u/dt;
h = dt;
timesteps_vec = 0:dt:time_u;

% transients are already done and following are the initial conditions
% after the transients.
ic = [-7.1778; -12.9972; 12.5330]; % initial conditions

u = zeros(max(size(ic)), nmax); % dynamics: u = [x; y; z]
u(:,1) = ic; % dynamics at first time-step

[m,nmax] = size(u); % m is the system dimension
M = m*nmax; 


% Gram-Schmidt parameters
nnorm1 = 1;
tnorm = nnorm1*dt;

cnorm = 10;
tcnorm = cnorm*dt;

wnorm = 1;
twnorm = wnorm*dt;

% c calculation parameters
nnorm2 = 10;

p = struct('sigma', sigma, 'r', r, 'b', b, 'nmax', nmax, 'dt', dt, 'm', m, ...
            'nnorm1', nnorm1, 'tnorm', tnorm, 'cnorm', cnorm, 'tcnorm', tcnorm, ...
            'wnorm', wnorm, 'twnorm', twnorm,'timesteps_vec', timesteps_vec);

tic
%% space-time  dynamics calculation

u = stdynamics(p,u,@rk4,@linLorenz63,@nonlinLorenz63);

% u = stdynamics(p,u,@etdrk4,@linLorenz63,@nonlinLorenz63);

% u = stdynamics(p,u,@etd1pc,@lorenz63rhs,@linLorenz63,@nonlinLorenz63);


%% dynamics plot
fig_attractor(u);

%% gs
% [v,R,laminst,lamgs] = gs(p,u,@rk4ts,@lor63jacobian);
% [v,R,laminst,lamgs] = gs(p,u,@rk4dyn,@lor63jacobian);
[v,R,laminst,lamgs] = gs(p,u,@etdrk4,@lor63jacobian);

%% laminst convergence fig
fig_laminstave(laminst);

%% functions
function u = stdynamics(p,u,solverfunc,linpartfunc,nonlinpartfunc)
% calculates the space time dynamics
linpart = linpartfunc(p);
nmax = p.nmax;
for n = 1:nmax
    u(:,n+1) = solverfunc(p,u(:,n),n,linpart,nonlinpartfunc);
end
end

function x = rk4(p,x,n,linpart,nonlinpartfunc)

% linpart = linpartfunc(p);
dynrhsfunc = @(p,x,linpart,nonlinpartfunc) linpart*x+nonlinpartfunc(x,n);
dt = p.dt;

% k values
k1 = dynrhsfunc(p, x , linpart, nonlinpartfunc);
k2 = dynrhsfunc(p, x + (0.5*dt)*k1 , linpart, nonlinpartfunc);
k3 = dynrhsfunc(p, x + (0.5*dt)*k2 , linpart, nonlinpartfunc);
k4 = dynrhsfunc(p, x + dt*k3 , linpart, nonlinpartfunc);
% dynamics
x = x + (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);

end


function u = etdrk4(p,u,n,linpart,nonlinpartfunc)
h = p.dt;
m = p.m;
t = p.timesteps_vec;

% precomputing some
L = linpart;
N = nonlinpartfunc;
invL = inv(L);

E1 = expm(L*h);
E2 = expm(L*h/2);
omega = invL * (E2-eye(m));

eta = h^(-2) * L^(-3);
alpha = eta * (- 4*eye(m) - L*h + E1 * (4*eye(m) - 3*L*h + (L*h)^2));
beta = eta * (2*eye(m) + L*h + E1 * (-2*eye(m) + L*h));
gamma = eta * (-4*eye(m) - 3*L*h - (L*h)^2 + E1 * (4*eye(m) - L*h));

%etd
an = E2 * u + omega * N(u,t(n));
bn = E2 * u + omega * N(an,t(n)+h/2);
cn = E2 * an + omega * ( 2 * N(bn,t(n)+h/2) - N(u,t(n)) );
u =  E1 * u + alpha * N(u, t(n)) ...
    + beta * 2 * ( N(an, t(n)+h/2) + N(bn, t(n)+h/2) ) ...
    + gamma * N(cn, t(n)+h) ;
end

function unp1 = etd1pc(p,un,n,linpart,nonlinpart)
L = linpart(p);
N = nonlinpart;
h = p.dt;
% predictor  step
un = diag(un);
uguess = un;
for iter = 2
    upred = expm(L*h)*un + 0.5*h*(N(uguess,n)+expm(L*h)*N(un,n));
    uguess = upred;
end
% corrector step
unp1 = expm(L*h)*un + 0.5*h*(N(upred,n)+expm(L*h)*N(un,n));
unp1 = diag(unp1);
end

function vec = lorenz63rhs(p,X,linpart,nonlinpart)
    L = linpart(p);
    vec = L*X + nonlinpart(X);
end

function L = linLorenz63(p)
sig = p.sigma;
r = p.r;
b = p.b;

L = [-sig, sig, 0;
        r, -1, 0;
        0, 0, -b];
end

function vec = nonlinLorenz63(un,n)
vec = [0; -un(1)*un(3); un(1)*un(2)];
end

function mat = lor63jacobian(p,Xt)
    sigma = p.sigma;
    r = p.r;
    b = p.b;
    mat = [-sigma sigma 0; ...
        (r-Xt(3)) -1 -Xt(1); ...
        Xt(2) Xt(1) -b];
end

function pertVecs = rk4ts(p,Xt,vt,jacobianFunc)
% pertVecs = rk4ts(p,Xt,vt,jacobianFunc)

dt = p.dt;
rhs = @(p,Xt,vt) jacobianFunc(p,Xt)*vt;

k1 = rhs(p,Xt,vt);
k2 = rhs(p, Xt, vt + dt*k1./2);
k3 = rhs(p, Xt, vt + dt*k2./2);
k4 = rhs(p, Xt, vt + dt*k3);

pertVecs = vt + (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);

end

function [v,R,laminst,lamgs] = gs(p,dynamics,solverfunc,jacobianFunc)
% [v,R,laminst,lamgs] = gs(p,dynamics,evolFunc,jacobianFunc)
% feed one time-step dynamics to the jacobian

nnorm1 = p.nnorm1;
tnorm = p.tnorm;

[M,nmax] = size(dynamics);
v(:,:,1) = eye(M); % initial perturbation vectors
R(:,:,1) = zeros(M,M);
laminst = zeros(M,1);
lamgs = zeros(M,1);

nonlinpart = @(v,n) [0; 0; 0];

fprintf("Calculating GS----------\n")
for n = 1:nmax
    % pertvecs = evolFunc(p,dynamics(:,n),v(:,:,n),jacobianFunc);
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