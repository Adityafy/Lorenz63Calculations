function u = etdCoxMatthewsRK4(p,u,Lfunc,Nfunc)
% u = etdCoxMatthewsRK4(p,u)

h = p.dt;
m = p.m;
nmax = p.nmax;

% precomputing some
L = Lfunc(p);
N = Nfunc;

E1 = expm(L*h);
E2 = expm(L*h/2);
% omega = inv(L) * (E2-eye(m));
omega = L \ (E2-eye(m));

eta = h^(-2) * L^(-3);
alpha = eta * (- 4 - L*h + E1 * (4 - 3*L*h + (L*h)^2));
beta = eta * 2 * (2 + L*h + E1 * (-2+L*h));
gamma = eta * (-4 - 3*L*h - (L*h)^2 + E1 * (4 - L*h));

% time integration
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

% time integration
for n = 1:nmax
    an = E2 * u(:,n) + omega * N(u(:,n));
    bn = E2 * u(:,n) + omega * N(an);
    cn = E2 * an + omega * ( 2 * N(bn) - N(u(:,n)) );

    u(:,n+1) =  E1 * u(:,n) + alpha * N(u(:,n)) ...
            + beta * ( N(an) + N(bn) ) ...
            + gamma * N(cn) ;

end

end