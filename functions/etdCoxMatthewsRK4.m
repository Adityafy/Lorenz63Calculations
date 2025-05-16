function u = etdCoxMatthewsRK4(p,u,dynrhsfunc,Lfunc,Nfunc)
% u = etdCoxMatthewsRK4(p,u,Lfunc,Nfunc)

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