%% Lyapunov exponents for Lorenz system
% Parameters and Variables

close all;
sigma = 10; %parameters
r = 28;
b = 8/3;

nmax = 100000; %number of iterations
deltat = 0.001; %time-step length
nnorm = 500;
tnorm = nnorm*deltat;

ic = [-1.6737 -2.6235 14.2779]; % initial conditions
x = [ic(1) zeros(1,nmax)]'; % trajectories
y = [ic(2) zeros(1,nmax)]';
z = [ic(3) zeros(1,nmax)]';

% Lorenz solver

j = 1;
while j<nmax+2
    x(j+1) = x(j) + deltat*(sigma*(y(j)-x(j)));
    y(j+1) = y(j) + deltat*(r*x(j)-y(j)-x(j)*z(j));
    z(j+1) = z(j) + deltat*(x(j)*y(j)-b*z(j));
    j = j+1;
end

figure(); % Plotting the Lorenz attractor
plot3(x(:,1),y(:,1),z(:,1), 'k',LineWidth=0.7);
xlabel('x'); ylabel('y'); zlabel('z'); view([0 0]);
%% Leading perturbation vector (with normalization)


ic1 = [0 0 1]; %initial conditions
dx1 = [ic1(1) zeros(1,nmax)]'; % trajectories
dy1 = [ic1(2) zeros(1,nmax)]';
dz1 = [ic1(3) zeros(1,nmax)]';

mag1 = zeros(1,nmax)'; %magnitude of the leading perturbation vector

lam1inst = [];
j = 1;
while j<=nmax+1
    dx1(j+1) = dx1(j) + deltat*(sigma*(dy1(j)-dx1(j)));
    dy1(j+1) = dy1(j) + deltat*(r*dx1(j)-dy1(j)-x(j)*dz1(j)-dx1(j)*z(j));
    dz1(j+1) = dz1(j) + deltat*(x(j)*dy1(j)+dx1(j)*y(j)-b*dz1(j));
    mag1(j+1) = sqrt(dx1(j+1)^2+dy1(j+1)^2+dz1(j+1)^2);
    if rem(j,nnorm) == 0 %renormalization
        lam1inst = [lam1inst (1/tnorm)*log(mag1(j+1))];
        norm1 = [dx1(j+1) dy1(j+1) dz1(j+1)]./norm([dx1(j+1) dy1(j+1) dz1(j+1)]);
        dx1(j+1) = norm1(1);
        dy1(j+1) = norm1(2);
        dz1(j+1) = norm1(3);
    end
    j = j+1;
end
d1 = [dx1 dy1 dz1];

figure(); %plotting lambda_1 inst wrt time
plot(lam1inst);
ylabel('\lambda_{1i} '); xlabel("time")
figure();
plot(mag1);

figure();
plot(cumsum(lam1inst)./(1:length(lam1inst)));
ylabel('Finite-time \lambda_1'); xlabel('time');
fprintf('Finite time 1st Lyapunov exponent is =')
disp(sum(lam1inst)/(nmax/nnorm));
%% 2nd perturbation vector

ic2 = [0 1 0]; %initial conditions
dx2 = [ic2(1) zeros(1,nmax)]'; % trajectories
dy2 = [ic2(2) zeros(1,nmax)]';
dz2 = [ic2(3) zeros(1,nmax)]';

mag2 = zeros(1,nmax)'; %magnitude of the leading perturbation vector

lam2inst = [];
j = 1;
while j<=nmax+1
    dx2(j+1) = dx2(j) + deltat*(sigma*(dy2(j)-dx2(j)));
    dy2(j+1) = dy2(j) + deltat*(r*dx2(j)-dy2(j)-x(j)*dz2(j)-dx2(j)*z(j));
    dz2(j+1) = dz2(j) + deltat*(x(j)*dy2(j)+dx2(j)*y(j)-b*dz2(j));
    mag2(j+1) = sqrt(dx2(j+1)^2+dy2(j+1)^2+dz2(j+1)^2);
    if rem(j,nnorm) == 0 %renormalization        
        d2star = [dx2(j+1) dy2(j+1) dz2(j+1)]-dot([dx2(j+1) dy2(j+1) dz2(j+1)], d1(j+1,:)).*d1(j+1,:);
        lam2inst = [lam2inst (1/tnorm)*log(norm(d2star))];
        norm2 = d2star./norm(d2star);
        dx2(j+1) = norm2(1);
        dy2(j+1) = norm2(2);
        dz2(j+1) = norm2(3);
    end
    j = j+1;
end
d2 = [dx2 dy2 dz2];

figure(); %plotting lambda_2 inst wrt time
plot(lam2inst);
ylabel('\lambda_{2i}'); xlabel("time")
figure();
plot(cumsum(lam2inst)./(1:length(lam2inst)));
ylabel('Finite-time \lambda_2'); xlabel('time');
fprintf('Finite time 2nd Lyapunov exponent is')
disp(sum(lam2inst)/(nmax/nnorm));
%% 3rd perturbation vector

ic3 = [1 0 0]; %initial conditions
dx3 = [ic3(1) zeros(1,nmax)]'; % trajectories
dy3 = [ic3(2) zeros(1,nmax)]';
dz3 = [ic3(3) zeros(1,nmax)]';

mag3 = zeros(1,nmax)'; %magnitude of the leading perturbation vector

lam3inst = [];
j = 1;
while j<=nmax+1
    dx3(j+1) = dx3(j) + deltat*(sigma*(dy3(j)-dx3(j)));
    dy3(j+1) = dy3(j) + deltat*(r*dx3(j)-dy3(j)-x(j)*dz3(j)-dx3(j)*z(j));
    dz3(j+1) = dz3(j) + deltat*(x(j)*dy3(j)+dx3(j)*y(j)-b*dz3(j));
    mag3(j+1) = sqrt(dx3(j+1)^2+dy3(j+1)^2+dz3(j+1)^2);
    if rem(j,nnorm) == 0 %renormalization        
        d3star = [dx3(j+1) dy3(j+1) dz3(j+1)]-...
            dot([dx3(j+1) dy3(j+1) dz3(j+1)], d2(j+1,:)).*d2(j+1,:)-...
            dot([dx3(j+1) dy3(j+1) dz3(j+1)], d1(j+1,:)).*d1(j+1,:);
        lam3inst = [lam3inst (1/tnorm)*log(norm(d3star))];
        norm3 = d3star./norm(d3star);
        dx3(j+1) = norm3(1);
        dy3(j+1) = norm3(2);
        dz3(j+1) = norm3(3);
    end
    j = j+1;
end

figure(); %plotting lambda_3 inst wrt time
plot(lam3inst);
ylabel('\lambda_{3i}'); xlabel("time")
figure();
plot(cumsum(lam3inst)./(1:length(lam3inst)));
ylabel('Finite-time \lambda_3'); xlabel('time');
fprintf('Finite time 3rd Lyapunov exponent is')
disp(sum(lam3inst)/(nmax/nnorm));