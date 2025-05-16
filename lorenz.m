% Lorenz CLV calculations runscript
close all;
clear all;
addpath('/Users/araj/Documents/Lorenz63Calculations/functions/');

%% parameters
sigma = 10; r = 28; b = 8/3;
%lorPara = [sigma, r, b]; % lorenz parameters

time_u = 200;

dt = 0.001;
nmax = time_u/dt;
h = dt;
timesteps_vec = 0:dt:time_u;

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

% u = stdynamics(p,u,@rk4dyn,@lorenz63rhs,@linLorenz63,@nonlinLorenz63);

u = stdynamics(p,u,@etdCoxMatthewsRK4,@lorenz63rhs,@linLorenz63,@nonlinLorenz63);


%% dynamics plot
fig_attractor(u);

%% gs
[v,R,laminst,lamgs] = gs(p,u,@rk4ts,@lor63jacobian);

%% laminst convergence fig
fig_laminstave(laminst);

%% c matrix
cmat = cmatrix(p,u,R,0,0);

%% CLVs
[wf,cle] = clv(p,u,v,cmat,@rk4ts,@lor63jacobian);

%%
fig_laminstave(cle);
title('CLEs');