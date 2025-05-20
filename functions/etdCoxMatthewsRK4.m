function u = etdCoxMatthewsRK4(p,u,n,linpart,nonlinpartfunc)
h = p.dt;
m = p.m;
t = p.timesteps_vec;

% precomputing some
L = linpart;
N = nonlinpartfunc;
invL = inv(L);

E1 = expm(L*h);
E2 = expm(L*h/2);
omega = invL * (E2-eye(m));

eta = h^(-2) * L^(-3);
alpha = eta * (- 4*eye(m) - L*h + E1 * (4*eye(m) - 3*L*h + (L*h)^2));
beta = eta * (2*eye(m) + L*h + E1 * (-2*eye(m) + L*h));
gamma = eta * (-4*eye(m) - 3*L*h - (L*h)^2 + E1 * (4*eye(m) - L*h));

%etd
an = E2 * u + omega * N(u,t(n));
bn = E2 * u + omega * N(an,t(n)+h/2);
cn = E2 * an + omega * ( 2 * N(bn,t(n)+h/2) - N(u,t(n)) );
u =  E1 * u + alpha * N(u, t(n)) ...
    + beta * 2 * ( N(an, t(n)+h/2) + N(bn, t(n)+h/2) ) ...
    + gamma * N(cn, t(n)+h) ;
end