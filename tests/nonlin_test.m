%% du/dt = u - u^3
% f4 : function 4 : rhs of du/dt = u - u^3

close all; clear all;
addpath('../functions/');

dt = [0.01 0.005 0.0025 0.001 0.0005];
% dt = 0.01;

for i = 1:length(dt)
    m = 1;
    time_u = 10;
    nmax = time_u/dt(i);
    timesteps_vec = linspace(0,time_u,nmax+1);

    p = struct('dt', dt(i), 'm', m, 'nmax', nmax, 'time_u', time_u, ...
        'timesteps_vec', timesteps_vec);

    %% du/dt = u - u^3
    % f3 : function 3 : rhs of du/dt = u-u^3
    dynrhsf3 = @(p,u,t) u-u^3;
    
    Lf3 = @(p) 1;

    Nf3 = @(u,t) -u^3;

    u(1) = 0.5;
    K = log(u(1)^2/sqrt(u(1)^2-1));

    for t = 1:nmax+1
        uexactf3(t) = exp(timesteps_vec(t))/sqrt(exp(2*timesteps_vec(t))+((1-u(1)^2)/u(1)^2));
    end
    uetdf3 = etdrk4(p,u,Lf3,Nf3);
    urk4f3 = rk4(p,u,dynrhsf3);
    ufoef3 = foe(p,u,dynrhsf3);
    uetd1f3 = etd1(p,u,Lf3,Nf3);
    uetd1expf3 = etd1explicit(p,u,Lf3,Nf3);

    etderrorf3(i) = mean(abs((uexactf3-uetdf3)./uexactf3));
    foeerrorf3(i) = mean(abs((uexactf3-ufoef3)./uexactf3));
    rk4errorf3(i) = mean(abs((uexactf3-urk4f3)./uexactf3));
    etd1errorf3(i) = mean(abs((uexactf3-uetd1f3)./uexactf3));
    etd1experrorf3(i) = mean(abs((uexactf3-uetd1expf3)./uexactf3));


end

%%
figure;
hold on;
plot(timesteps_vec,uetdf3,'-o','color','blue','DisplayName','u etdcm');
plot(timesteps_vec,urk4f3,'-s','color','red', 'DisplayName','u rk4');
plot(timesteps_vec,uexactf3,'color','black', 'DisplayName','u exact','LineWidth',2);
hold off;
title('du/dt = u - u^3');
legend;

%%
figure;
% plot(log10(dt),log10(etderrorf3),'-o','DisplayName','ETDRK4','LineWidth',1);
hold on;
plot(log10(dt),log10(foeerrorf3), '-s','color', 'black','DisplayName','First order Euler', ...
    'LineWidth',1);
plot(log10(dt),log10(rk4errorf3), '-^','color', 'red','DisplayName', ...
    'RK4','LineWidth',1);
% plot(log10(dt),log10(etd1errorf3), '-^','color', 'green','DisplayName', ...
%     'ETD1','LineWidth',1);
plot(log10(dt),log10(etd1experrorf3), '-v','color', 'blue','DisplayName', ...
    'ETD1exp','LineWidth',1);
hold off;
legend;
xlim([log10(min(dt)) log10(max(dt))]);
subtitle('Error convergence for du/dt = u - u^3');
set(gca,'TickLabelInterpreter','tex','FontSize',15, ...
    'XMinorTick','on','YMinorTick','on','Box','on');
axis square;
xlabel('$\log{(\Delta t)}$','Interpreter','latex','FontSize',20);
ylabel('$\log{(E)}$','Interpreter','latex','FontSize',20);
subtitle('Error convergence for du/dt = u - u^3');




%% FOE function
function u = foe(p,u,rhs)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
dt = p.dt;
nmax = p.nmax;
t = p.timesteps_vec;
for n = 1:nmax
    k1 = rhs(p, u(:,n), t(n));
    u(:,n+1) = u(:,n) + dt * k1;
end
end


%% RK4 function
function u = rk4(p,u,rhs)

h = p.dt;
nmax = p.nmax;
t = p.timesteps_vec;

for n = 1:nmax
    % k values
    k1 = rhs(p, u(:,n), t(n));
    k2 = rhs(p, u(:,n) + (0.5*h)*k1, t(n)+h/2);
    k3 = rhs(p, u(:,n) + (0.5*h)*k2, t(n)+h/2);
    k4 = rhs(p, u(:,n) + h*k3, t(n)+h);
    % dynamics
    u(:,n+1) = u(:,n) + (1/6) * h * (k1 + 2*k2 + 2*k3 + k4);
end
end


%% ETD1 function
function u = etd1(p,u,Lfunc,Nfunc)

dt = p.dt;
nmax = p.nmax;
t = p.timesteps_vec;
L = Lfunc(p);
N = Nfunc;

for n = 1:nmax
    u(:,n+1) = u(:,n)*exp(L*dt) + N(u(:,n),t(n)) * (exp(L*dt)-1)/L;
end

end

%% ETD1 function
function u = etd1explicit(p,u,Lfunc,Nfunc)
% This is done by the predictor corrector method used in Cross et al.
% (1994). There is an extra loop in the time integration loop so that
% predictor guess is converged.
dt = p.dt;
nmax = p.nmax;
t = p.timesteps_vec;
L = Lfunc(p);
N = Nfunc;

for n = 1:nmax
    % predictor1
    u_pred = 0.25; % initial guess for predictor step
    for i = 1:200
        u_pred = u(:,n)*expm(L*dt) + 0.5*dt*N(u_pred,t(n)) + ...
                        0.5*dt*expm(L*dt)*N(u(:,n),t(n));
    end
    % corrector (final step for u_{n+1})
    u(:,n+1) = u(:,n)*expm(L*dt) + 0.5*dt*N(u_pred,t(n)) + ...
                        0.5*dt*expm(L*dt)*N(u(:,n),t(n));
end

end

%% ETD RK4 function
% function u = etdrk4(p,u,Lfunc,Nfunc)
% % u = etdCoxMatthewsRK4(p,u)
% 
% h = p.dt;
% m = p.m;
% nmax = p.nmax;
% t = p.timesteps_vec;
% 
% % precomputing some
% L = Lfunc(p);
% N = Nfunc;
% 
% % E1 = expm(L*h);
% % E2 = expm(L*h/2);
% E1 = exp(L*h);
% E2 = exp(L*h/2);
% % omega = inv(L) * (E2-eye(m));
% omega = L \ (E2-eye(m));
% 
% eta = h^(-2) * L^(-3);
% alpha = eta * (- 4 - L*h + E1 * (4 - 3*L*h + (L*h)^2));
% beta = eta * (2 + L*h + E1 * (-2+L*h));
% gamma = eta * (-4 - 3*L*h - (L*h)^2 + E1 * (4 - L*h));
% 
% % time integration
% for n = 1:nmax
%     an = E2 * u(:,n) + omega * N(u(:,n),t(n));
%     bn = E2 * u(:,n) + omega * N(an,t(n)+h/2);
%     cn = E2 * an + omega * ( 2 * N(bn,t(n)+h/2) - N(u(:,n),t(n)) );
% 
%     u(:,n+1) =  E1 * u(:,n) + alpha * N(u(:,n), t(n)) ...
%         + beta * 2 * ( N(an, t(n)+h/2) + N(bn, t(n)+h/2) ) ...
%         + gamma * N(cn, t(n)+h) ;
% 
% end
% 
% end

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
