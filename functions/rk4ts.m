function pertVecs = rk4ts(para,h,Xt,vt,jacobianFunc)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)


k1 = jacobianFunc(para,Xt,vt);
% k2 = jacobianFunc(Xt, d1(:,t)+h*k1(:,1)/2, d2(:,t)+h*k1(:,2)/2, d3(:,t)+h*k1(:,3)/2);
k2 = jacobianFunc(para, Xt, vt + h*k1./2);
% k3 = jacobianFunc(x(t),y(t),z(t), d1(:,t)+h*k2(:,1)/2, d2(:,t)+h*k2(:,2)/2, d3(:,t)+h*k2(:,3)/2);
k3 = jacobianFunc(para, Xt, vt + h*k2./2);
k4 = jacobianFunc(para, Xt, vt + h*k3);

% rk4 applied to perturbation vectors
% A = [d1(:,t) d2(:,t) d3(:,t)] + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
pertVecs = vt + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);

end
