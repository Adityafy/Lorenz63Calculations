% rossler

close all;
addpath('/Users/araj/Documents/Lorenz63Calculations/functions/');

%% parameters
a = 0.15; b = 0.20; c = 10;


trans_time_u = 501;
time_u = 1000;

dt = 0.001;
nmax_trans = trans_time_u/dt;
nmax = time_u/dt;
h = dt;

ic = [0.1; 0.1; 0.1]; % initial conditions
trans_dynamics = zeros(max(size(ic)), nmax_trans);
trans_dynamics(:,1) = ic;
dynamics = zeros(max(size(ic)), nmax);


[m,nmax] = size(dynamics); % m is the system dimension
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

p = struct('a', a, 'b', b, 'c', c, 'nmax', nmax, 'dt', dt, 'm', m, ...
            'nnorm1', nnorm1, 'tnorm', tnorm, 'cnorm', cnorm, 'tcnorm', tcnorm, ...
            'wnorm', wnorm, 'twnorm', twnorm);

tic
%% dynamics calculation
trans_dynamics = rk4dyn(p,dynamics,@rosslerrhs);
dynamics(:,1) = trans_dynamics(:,end);
clear trans_dynamics;
dynamics = rk4dyn(p,dynamics,@rosslerrhs);

%% dynamics plot
fig_attractor(dynamics);

%% gs
[v,R,laminst,lamgs] = gs(p,dynamics,@rk4ts,@rosslerjacobian);

%% laminst convergence fig
fig_laminstave(laminst);

%% c matrix
cmat = cmatrix(p,dynamics,R,0,0);

%% CLVs
[wf,cle] = clv(p,dynamics,v,cmat,@rk4ts,@rosslerjacobian);

