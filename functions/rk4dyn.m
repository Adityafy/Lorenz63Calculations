function dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
% dynamics = rk4dyn(para,dynamics,dynFunc,h,total_time)
for j = 1:total_time
    % k values
    k1 = dynFunc(para, dynamics(:,j));
    k2 = dynFunc(para, dynamics(:,j) + (0.5*h)*k1);
    k3 = dynFunc(para, dynamics(:,j) + (0.5*h)*k2);
    k4 = dynFunc(para, dynamics(:,j) + h*k3);
    % dynamics
    dynamics(:,j+1) = dynamics(:,j) + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
end
end

