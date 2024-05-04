% Lorenz CLV calculations runscript
close all;
clear all;
addpath('/Users/araj/Documents/Lorenz63Calculations/functions/');

%% parameters
sigma = 10; r = 28; b = 8/3;
%lorPara = [sigma, r, b]; % lorenz parameters

time_u = 500;

dt = 0.001;
nmax = time_u/dt;
h = dt;

ic = [-7.1778; -12.9972; 12.5330]; % initial conditions
dynamics = zeros(max(size(ic)), nmax);
dynamics(:,1) = ic;

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

p = struct('sigma', sigma, 'r', r, 'b', b, 'nmax', nmax, 'dt', dt, 'm', m, ...
            'nnorm1', nnorm1, 'tnorm', tnorm, 'cnorm', cnorm, 'tcnorm', tcnorm, ...
            'wnorm', wnorm, 'twnorm', twnorm);

tic
%% dynamics calculation
dynamics = rk4dyn(p,dynamics,@lorenz63rhs);

%% dynamics plot
fig_attractor(dynamics);

%% gs
[v,R,laminst,lamgs] = gs(p,dynamics,@rk4ts,@lor63jacobian);

%% laminst convergence fig
fig_laminstave(laminst);

%% c matrix
cmat = cmatrix(p,dynamics,R,0,0);

%% CLVs
[wf,cle] = clv(p,dynamics,v,cmat,@rk4ts,@lor63jacobian);

%%
figure();
hold on;
%timeaxix = linspace(createw*deltat,nmax*deltat,deltat*(nmax-createw)/nnorm3);
plot(cumsum(cle(1,:))./(1:length(cle(1,:))), 'LineWidth',1);
set(gca,'TickLabelInterpreter','latex');
ylabel('Finite-time CLE $\lambda_1$ , $\lambda_2$ , $\lambda_3$', 'Interpreter','latex', 'FontWeight','bold');
xlabel('Time-steps','Interpreter','latex', 'FontWeight','bold');
plot(cumsum(cle(2,:))./(1:length(cle(2,:))), 'LineWidth',1);
plot(cumsum(cle(3,:))./(1:length(cle(3,:))), 'LineWidth',1);
legend('\lambda_1','\lambda_2','\lambda_3');
% xlim([0 length(cle(1,:))]);
% xticks([createw*deltat deltat*(nmax-createw)/length(cle(1,:)) nmax*deltat]);
hold off;