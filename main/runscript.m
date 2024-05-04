
clear all;
addpath('../functions/');

%% parameters
sigma = 10; r = 28; b = 8/3;
%lorPara = [sigma, r, b]; % lorenz parameters

time_u = 200;

dt = 0.01;
nmax = time_u/dt;
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
            'nnorm1', nnorm1, 'nnorm2', nnorm2, 'tnorm', tnorm);

%% exponential time integration (Cox and Matthews)
u = etdCoxMatthewsRK4(p,u,@linLorenz63,@nonlinLorenz63);


%% dynamics plot
fig_attractor(u);