close all;
addpath('/Users/araj/Documents/Lorenz63Calculations/functions/')
sigma = 10; r = 28; b = 8/3;
lorPara = [sigma, r, b]; % lorenz parameters

total_time = 200000;
deltat = 0.001;
h = deltat;

ic = [-7.1778; -12.9972; 12.5330]; % initial conditions
dynamics = ic;

[m,n] = size(dynamics);
M = m*n; % system dimension

% dynamics calculation
dynamics = rk4dyn(lorPara,dynamics,@lorenz63,h,total_time);

% dynamics plot
fig_attractor(dynamics);

% Gram-Schmidt
nnorm = 100;
tnorm = nnorm*deltat;
[v,R,laminst,lamgs] = gs(lorPara,h,nnorm,dynamics,@rk4ts,@lor63jacobian);

% laminst convergence fig
fig_laminstave(laminst);

