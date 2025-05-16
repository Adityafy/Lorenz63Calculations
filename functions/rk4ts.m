function pertVecs = rk4ts(p,Xt,vt,jacobianFunc)
% pertVecs = rk4ts(p,Xt,vt,jacobianFunc)

dt = p.dt;
rhs = @(p,Xt,vt) jacobianFunc(p,Xt)*vt;

k1 = rhs(p,Xt,vt);
k2 = rhs(p, Xt, vt + dt*k1./2);
k3 = rhs(p, Xt, vt + dt*k2./2);
k4 = rhs(p, Xt, vt + dt*k3);

pertVecs = vt + (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);

end