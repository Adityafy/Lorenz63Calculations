% testing error convergence for the numerical solution of
% du/dt = u and dv/dt = - u^3
% using etd (cox and matthews), rk4, first order Euler

close all; clear all;
addpath('../functions/');

dt = [0.5 0.2 0.1 0.05 0.01 0.005 0.001 0.0005];
% dt = 10.^(0:-1:-4);
% dt = 1;
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

    uetdf1 = etdrk4(p,u,Lf1,Nf1);
    urk4f1 = rk4(p,u,dynrhsf1);
    uexactf1 = uexactf1fn(p,u);
    ufoef1 = foedyn(p,u,dynrhsf1);

    % errors
    % etderrorf1(i) = mean(abs(uexactf1-uetdf1(1:end-1)));
    % foeerrorf1(i) = mean(abs(uexactf1-ufoef1(1:end-1)));
    % rk4errorf1(i) = mean(abs(uexactf1-urk4f1(1:end-1)));

    % etderrorf1(i) = mean(abs(uexactf1-uetdf1));
    % foeerrorf1(i) = mean(abs(uexactf1-ufoef1));
    % rk4errorf1(i) = mean(abs(uexactf1-urk4f1));

    etderrorf1(i) = mean(abs((uexactf1-uetdf1)./uexactf1));
    foeerrorf1(i) = mean(abs((uexactf1-ufoef1)./uexactf1));
    rk4errorf1(i) = mean(abs((uexactf1-urk4f1)./uexactf1));

    %% du/dt = -u^3
    % f2 : function 2 : rhs of du/dt = -u^3
    dynrhsf2 = @(p,u) -u^3;
    uexactf2fn = @(p,u) 1./sqrt(2*p.time+(1/u(1)^2));

    Lf2 = @(p) 0;

    Nf2 = @(u) -u^3;

    % u = zeros(1,nmax);
    u(1) = 0.5;

    uetdf2 = etdrk4(p,u,Lf2,Nf2);
    urk4f2 = rk4(p,u,dynrhsf2);
    uexactf2 = uexactf2fn(p,u);
    ufoef2 = foedyn(p,u,dynrhsf2);

    % errors
    % etderrorf2(i) = mean(abs(uexactf2-uetdf2(1:end-1)));
    % foeerrorf2(i) = mean(abs(uexactf2-ufoef2(1:end-1)));
    % rk4errorf2(i) = mean(abs(uexactf2-urk4f2(1:end-1)));

    % etderrorf2(i) = mean(abs(uexactf2-uetdf2));
    % foeerrorf2(i) = mean(abs(uexactf2-ufoef2));
    % rk4errorf2(i) = mean(abs(uexactf2-urk4f2));

    etderrorf2(i) = mean(abs((uexactf2-uetdf2)./uexactf2));
    foeerrorf2(i) = mean(abs((uexactf2-ufoef2)./uexactf2));
    rk4errorf2(i) = mean(abs((uexactf2-urk4f2)./uexactf2));

    %% du/dt = u - u^3
    % f3 : function 3 : rhs of du/dt = u-u^3
    dynrhsf3 = @(p,u) u-u^3;
    
    Lf3 = @(p) 1;

    Nf3 = @(u) -u^3;

    % u = zeros(1,nmax);
    u(1) = 0.5;
    K = log(u(1)^2/sqrt(u(1)^2-1));

    for t = 1:nmax+1
        uexactf3(t) = exp(time(t))/sqrt(exp(2*time(t))+((1-u(1)^2)/u(1)^2));
    end
    uetdf3 = etdrk4(p,u,Lf3,Nf3);
    urk4f3 = rk4(p,u,dynrhsf3);
    % uexactf3 = uexactf3fn(K,p,u);
    ufoef3 = foedyn(p,u,dynrhsf3);

    % errors
    % etderrorf3(i) = mean(abs(uexactf3-uetdf3(1:end-1)));
    % foeerrorf3(i) = mean(abs(uexactf3-ufoef3(1:end-1)));
    % rk4errorf3(i) = mean(abs(uexactf3-urk4f3(1:end-1)));

    % etderrorf3(i) = mean(abs(uexactf3-uetdf3));
    % foeerrorf3(i) = mean(abs(uexactf3-ufoef3));
    % rk4errorf3(i) = mean(abs(uexactf3-urk4f3));

    etderrorf3(i) = mean(abs((uexactf3-uetdf3)./uexactf3));
    foeerrorf3(i) = mean(abs((uexactf3-ufoef3)./uexactf3));
    rk4errorf3(i) = mean(abs((uexactf3-urk4f3)./uexactf3));


        %% du/dt = cu + sin(t)
    % f4 : function 4 : rhs of du/dt = cu + sin(t)
    c = -100;
    dynrhsf4 = @(p,u) c*u + sin(p.time);
    
    Lf4 = @(p) c;

    Nf4 = @(p) sin(p.time);

    % u = zeros(1,nmax);
    u(1) = 0.2;
    

    for t = 1:nmax+1
        uexactf4(t) = u(1)*exp(c*t)+(exp(c*t)-c*sin(t)-cos(t))/(1+c^2);
    end
    uetdf4 = etdrk4(p,u,Lf4,Nf4);
    urk4f4 = rk4(p,u,dynrhsf4);
    % uexactf3 = uexactf3fn(K,p,u);
    ufoef4 = foedyn(p,u,dynrhsf4);

    % errors
    % etderrorf3(i) = mean(abs(uexactf3-uetdf3(1:end-1)));
    % foeerrorf3(i) = mean(abs(uexactf3-ufoef3(1:end-1)));
    % rk4errorf3(i) = mean(abs(uexactf3-urk4f3(1:end-1)));

    % etderrorf3(i) = mean(abs(uexactf3-uetdf3));
    % foeerrorf3(i) = mean(abs(uexactf3-ufoef3));
    % rk4errorf3(i) = mean(abs(uexactf3-urk4f3));

    etderrorf4(i) = mean(abs((uexactf4-uetdf4)./uexactf4));
    foeerrorf4(i) = mean(abs((uexactf4-ufoef4)./uexactf4));
    rk4errorf4(i) = mean(abs((uexactf4-urk4f4)./uexactf4));


end
%%
figure;
% subplot(2,2,1)
hold on;
% plot(time,ufoef1(1:end-1),'o', 'DisplayName','u, First order Euler');
plot(time,uetdf1,'-o','color','blue','DisplayName','u etdcm');
plot(time,urk4f1,'color','red', 'DisplayName','u rk4','LineWidth',2);
plot(time,uexactf1,'-', 'DisplayName','u, Exact', 'LineWidth',1.5);
hold off;
title('du/dt = -u');
legend;

%%
figure;
% subplot(2,2,2);
plot(log10(dt),log10(etderrorf1),'-o','DisplayName','ETDRK4 error'); 
hold on;
plot(log10(dt),log10(foeerrorf1), '-s','DisplayName','First order Euler error');
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

%%
figure;
% subplot(2,2,3);
hold on;
plot(time,ufoef2,'o', 'DisplayName','u, First order Euler');
% plot(time,uetdf2(1:end-1),'-o','color','blue','DisplayName','u etdcm');
% plot(time,urk4f2(1:end-1),'color','red', 'DisplayName','u rk4');
plot(time,uexactf2,'-', 'DisplayName','u, Exact','LineWidth',1.5);
hold off;
title('du/dt = - u^3');
legend;

%%
figure;
% subplot(2,2,4);
% subplot(2,2,2);
plot(log10(dt),log10(etderrorf2),'--','DisplayName','ETDRK4 error'); 
hold on;
plot(log10(dt),log10(foeerrorf2), '-o','DisplayName','First order Euler error');
plot(log10(dt),log10(rk4errorf2), '-^','color', 'black','DisplayName','RK4 error');
hold off;
% ylim([-10 -1]);
% set(gca,'XDir','reverse');
legend;
% set(gca,'YScale','log');
% set(gca,'XScale','log');
xlabel('$\log{(dt)}$','Interpreter','latex');
ylabel('$\log{(E)}$','Interpreter','latex');
set(gca,'TickLabelInterpreter','tex','FontSize',15);
subtitle('Error convergence for du/dt = -u^3');

%%
figure;
% subplot(2,2,3);
hold on;
% plot(time,ufoef3,'o', 'DisplayName','u, First order Euler');
plot(time,uetdf3,'-o','color','blue','DisplayName','u etdcm');
plot(time,urk4f3,'-s','color','red', 'DisplayName','u rk4');
plot(time,uexactf3,'color','black', 'DisplayName','u exact','LineWidth',2);
% legend;
% plot(time,uexactf3,'-', 'DisplayName','u, Exact','LineWidth',1.5);
hold off;
title('du/dt = u - u^3');
legend;

%%
figure;
% subplot(2,2,4);
% subplot(2,2,2);
plot(log10(dt),log10(etderrorf3),'-o','DisplayName','ETDRK4'); 
hold on;
plot(log10(dt),log10(foeerrorf3), '-s','DisplayName','First order Euler');
plot(log10(dt),log10(rk4errorf3), '-^','color', 'black','DisplayName','RK4');
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
subtitle('Error convergence for du/dt = u - u^3');

% %%
% figure;
% hold on;
% plot(time,uetdf4,'-o','color','blue','DisplayName','u etdcm');
% plot(time,urk4f4,'-s','color','red', 'DisplayName','u rk4');
% plot(time,uexactf4,'color','black', 'DisplayName','u exact','LineWidth',2);
% hold off;
% title('du/dt = cu + sin(t)');
% legend;
% 
% %%
% figure;
% plot(log10(dt),log10(etderrorf4),'-o','DisplayName','ETDRK4'); 
% hold on;
% plot(log10(dt),log10(foeerrorf4), '-s','DisplayName','First order Euler');
% plot(log10(dt),log10(rk4errorf4), '-^','color', 'black','DisplayName','RK4');
% hold off;
% legend;
% xlabel('$\log{(\Delta t)}$','Interpreter','latex');
% ylabel('$\log{(E)}$','Interpreter','latex');
% subtitle('Error convergence for du/dt = cu + sin(t)');
% set(gca,'TickLabelInterpreter','tex','FontSize',15);
% subtitle('Error convergence for du/dt = cu + sin(t)');

%% RK4 function
function u = rk4(p,u,rhs)

dt = p.dt;
nmax = p.nmax;
for j = 1:nmax
    % k values
    k1 = rhs(p, u(:,j));
    k2 = rhs(p, u(:,j) + (0.5*dt)*k1);
    k3 = rhs(p, u(:,j) + (0.5*dt)*k2);
    k4 = rhs(p, u(:,j) + dt*k3);
    % dynamics
    u(:,j+1) = u(:,j) + (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);
end
end


%% ETD RK4 function
function u = etdrk4(p,u,Lfunc,Nfunc)
% u = etdCoxMatthewsRK4(p,u)

h = p.dt;
m = p.m;
nmax = p.nmax;

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
    an = E2 * u(:,n) + omega * N(u(:,n));
    bn = E2 * u(:,n) + omega * N(an);
    cn = E2 * an + omega * ( 2 * N(bn) - N(u(:,n)) );

    u(:,n+1) =  E1 * u(:,n) + alpha * N(u(:,n)) ...
            + beta * 2 * ( N(an) + N(bn) ) ...
            + gamma * N(cn) ;

end

% % time integration
% for n = 1:nmax
%     an = E2 * u(:,n) + omega * N(u(:,n));
%     bn = E2 * u(:,n) + omega * N(an+(h/2)*an);
%     cn = E2 * an + omega * ( 2 * N(bn+(h/2)*bn) - N(u(:,n)) );
% 
%     u(:,n+1) =  E1 * u(:,n) + alpha * N(u(:,n)) ...
%             + beta * ( N(an+(h/2)*an) + N(bn+(h/2)*bn) ) ...
%             + gamma * N(cn+h*cn) ;
% 
% end

% % time integration
% for n = 1:nmax
%     an = E2 * u(:,n) + omega * N(u(:,n));
%     bn = E2 * u(:,n) + omega * N(u(:,n)+(h/2)*an);
%     cn = E2 * an + omega * ( 2 * N(u(:,n)+(h/2)*bn) - N(u(:,n)) );
% 
%     u(:,n+1) =  E1 * u(:,n) + alpha * N(u(:,n)) ...
%             + beta * ( N(u(:,n)+(h/2)*an) + N(u(:,n)+(h/2)*bn) ) ...
%             + gamma * N(u(:,n)+h*cn) ;
% 
% end

end