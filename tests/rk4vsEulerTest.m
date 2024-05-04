% testing error convergence for the numerical solution of
% du/dt = u and dv/dt = - u^3
% using etd (cox and matthews), rk4, first order Euler

close all; clear all;
addpath('../functions/');

dt = [0.5 0.2 0.1 0.05 0.01 0.005 0.001 0.0005];
% dt = 10.^(0:-1:-5);
% dt = 10^-4;
for i = 1:length(dt)
    % dt = 0.01;
    m = 1;
    time_u = 10;
    nmax = time_u/dt(i);
    time = linspace(0,time_u,nmax+1);

    p = struct('dt', dt(i), 'm', m, 'nmax', nmax, 'time_u', time_u, ...
        'time', time);

    %% du/dt = -u
    % f1 : function 1: rhs of du/dt = -u
    dynrhsf1 = @(p,u) -u;
    uexactf1fn = @(p,u) u(1)*exp(-p.time);

    Lf1 = @(p) -1;

    Nf1 = @(u) 0;

    u = zeros(1,round(nmax));
    u(1) = 1;

    % uetdf1 = etdrk4(p,u,Lf1,Nf1);
    urk4f1 = rk4dyn(p,u,dynrhsf1);
    uexactf1 = uexactf1fn(p,u);
    ufoef1 = foedyn(p,u,dynrhsf1);

    % etderrorf1(i) = mean(abs(uexactf1-uetdf1));
    foeerrorf1(i) = mean(abs(uexactf1-ufoef1));
    rk4errorf1(i) = mean(abs(uexactf1-urk4f1));


end
%%
figure;
% subplot(2,2,1)
hold on;
plot(time,ufoef1,'o', 'DisplayName','u, First order Euler');
plot(time,urk4f1,'^','color','black', 'DisplayName','u rk4');
plot(time,uexactf1,'-', 'DisplayName','u, Exact', 'LineWidth',1.5);
hold off;
set(gca,'TickLabelInterpreter','tex','FontSize',15);
xlabel('$t$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
title('du/dt = -u');
legend;

%%
figure;
% subplot(2,2,2);
% plot(dt,etderrorf1,'--','DisplayName','ETDRK4 error'); 
hold on;
plot(log10(dt),log10(foeerrorf1), '-o','DisplayName','First order Euler error');
plot(log10(dt),log10(rk4errorf1), '-^','color', 'black','DisplayName','RK4 error');
hold off;
% ylim([-10 -1]);
% set(gca,'XDir','reverse');
legend;
% set(gca,'YScale','log');
% set(gca,'XScale','log');
xlabel('$\log{(dt)}$','Interpreter','latex');
ylabel('$\log{(E)}$','Interpreter','latex');
subtitle('Error convergence for du/dt = -u');
set(gca,'TickLabelInterpreter','tex','FontSize',15);



% %% RK4 function
% function u = rk4(p,u,rhs)
% 
% dt = p.dt;
% nmax = p.nmax;
% for j = 1:nmax
%     % k values
%     k1 = rhs(p, u(:,j));
%     k2 = rhs(p, u(:,j) + (0.5*dt)*k1);
%     k3 = rhs(p, u(:,j) + (0.5*dt)*k2);
%     k4 = rhs(p, u(:,j) + dt*k3);
%     % dynamics
%     u(:,j+1) = u(:,j) + (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);
% end
% end
% 
% %% FOE function
% function u = foe(p,u,rhs)
% dt = p.dt;
% nmax = p.nmax;
% for j = 1:nmax
%     k1 = rhs(p, u(:,j));
%     u(:,j+1) = u(:,j) + dt*k1;
% end
% end

