%% du/dt = cu + sin(t)
% f4 : function 4 : rhs of du/dt = cu + sin(t)

close all; clear all;
addpath('../functions/');

dt = [0.01 0.005 0.0025 0.001 0.0005];

for i = 1:length(dt)
    m = 1;
    time_u = 10;
    nmax = time_u/dt(i);
    time = linspace(0,time_u,nmax+1);

    p = struct('dt', dt(i), 'm', m, 'nmax', nmax, 'time_u', time_u, ...
        'time', time);

    c = -100;
    dynrhsf4 = @(p,u,t) c*u + sin(t);

    Lf4 = @(p) c;

    Nf4 = @(u,t) sin(t);

    % u = zeros(1,nmax);
    u(1) = 10;


    for t = 1:nmax+1
        uexactf4(t) = u(1)*exp(c*time(t)) + ...
            (exp(c*time(t))-c*sin(time(t))-cos(time(t)))/(1+c^2);
    end

    uetdf4 = etdrk4(p,u,Lf4,Nf4);
    urk4f4 = rk4(p,u,dynrhsf4);
    ufoef4 = foe(p,u,dynrhsf4);
    uetd1 = etd1(p,u,Lf4,Nf4);

    etderrorf4(i) = mean(abs((uexactf4-uetdf4)./uexactf4));
    foeerrorf4(i) = mean(abs((uexactf4-ufoef4)./uexactf4));
    rk4errorf4(i) = mean(abs((uexactf4-urk4f4)./uexactf4));
    etd1error(i) = mean(abs((uexactf4-uetd1)./uexactf4));
end

%%
figure;
hold on;
plot(time,uetdf4,'-o','color','blue','DisplayName','u etdcm');
plot(time,urk4f4,'-s','color','red', 'DisplayName','u rk4');
plot(time,uexactf4,'color','black', 'DisplayName','u exact','LineWidth',2);
hold off;
title('du/dt = cu + sin(t)');
legend;

%%
figure;
plot(log10(dt),log10(etderrorf4),'-o','DisplayName','ETDRK4');
hold on;
plot(log10(dt),log10(foeerrorf4), '-s','DisplayName','First order Euler');
plot(log10(dt),log10(rk4errorf4), '-^','color', 'black','DisplayName','RK4');
plot(log10(dt),log10(etd1error), '-^','color', 'green','DisplayName','ETD1');
hold off;
legend;
xlabel('$\log{(\Delta t)}$','Interpreter','latex');
ylabel('$\log{(E)}$','Interpreter','latex');
subtitle('Error convergence for du/dt = cu + sin(t)');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
subtitle('Error convergence for du/dt = cu + sin(t)');


%% Functions

% FOE function
function u = foe(p,u,rhs)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
dt = p.dt;
nmax = p.nmax;
t = p.time;
for n = 1:nmax
    k1 = rhs(p, u(:,n), t(n));
    u(:,n+1) = u(:,n) + dt * k1;
end
end


% RK4 function
function u = rk4(p,u,rhs)

dt = p.dt;
nmax = p.nmax;
t = p.time;

for n = 1:nmax
    % k values
    k1 = rhs(p, u(:,n), t(n));
    k2 = rhs(p, u(:,n) + (0.5*dt)*k1, t(n)+dt/2);
    k3 = rhs(p, u(:,n) + (0.5*dt)*k2, t(n)+dt/2);
    k4 = rhs(p, u(:,n) + dt*k3, t(n)+dt);
    % dynamics
    u(:,n+1) = u(:,n) + (1/6) * dt * (k1 + 2*k2 + 2*k3 + k4);
end
end

% ETD1 function
function u = etd1(p,u,Lfunc,Nfunc)

dt = p.dt;
nmax = p.nmax;
t = p.time;
L = Lfunc(p);
N = Nfunc;

for n = 1:nmax
    u(:,n+1) = u(:,n)*exp(L*dt) + N(u(:,n),t(n)) * (exp(L*dt)-1)/L;
end

end

% ETD RK4 function
function u = etdrk4(p,u,Lfunc,Nfunc)
% u = etdCoxMatthewsRK4(p,u)

h = p.dt;
m = p.m;
nmax = p.nmax;
t = p.time;

% precomputing some
L = Lfunc(p);
N = Nfunc;

% E1 = expm(L*h);
% E2 = expm(L*h/2);
E1 = exp(L*h);
E2 = exp(L*h/2);
% omega = inv(L) * (E2-eye(m));
omega = L \ (E2-eye(m));

eta = h^(-2) * L^(-3);
alpha = eta * (- 4 - L*h + E1 * (4 - 3*L*h + (L*h)^2));
beta = eta * (2 + L*h + E1 * (-2+L*h));
gamma = eta * (-4 - 3*L*h - (L*h)^2 + E1 * (4 - L*h));

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