%% Lorenz with RK and QR
close all;
% clear all;

%% Lorenz solver with QR
sigma = 10; %parameters
r = 28;
b = 8/3;

nmax = 200000; %number of iterations

nnorm = 1000;
deltat = 0.001; %time-step length
tnorm = nnorm*deltat;

% ic = [1; 1; 1;];
ic = [-7.1778; -12.9972; 12.5330]; % initial conditions
x = [ic(1) zeros(1,nmax)].'; % trajectories
y = [ic(2) zeros(1,nmax)].';
z = [ic(3) zeros(1,nmax)].';

% d1 = [dx1;
%       dy1;
%       dz1] (similar for d2)
ic1 = [0; 0; 1]; ic2 = [0; 1; 0]; ic3 = [1; 0; 0];
d1 = [ic1 zeros(3, nmax)];
d2 = [ic2 zeros(3, nmax)];
d3 = [ic3 zeros(3, nmax)];

h = deltat;
lam1inst = [];
lam2inst = [];
lam3inst = [];

% mag1 = [];
% mag2 = [];
% mag3 = [];

j = 1;
while j<nmax+2

    % k values, rk4 for lorenz
    k1 = lorenz_rhs(x(j), y(j), z(j));
    k2 = lorenz_rhs(x(j)+h*k1(1)/2, y(j)+h*k1(2)/2, z(j)+h*k1(3)/2);
    k3 = lorenz_rhs(x(j)+h*k2(1)/2, y(j)+h*k2(2)/2, z(j)+h*k2(3)/2);
    k4 = lorenz_rhs(x(j)+h*k3(1), y(j)+h*k3(2), z(j)+h*k3(3));

    % rk4 applied to lorenz system
    x(j+1) = x(j) + (1/6)*h*(k1(1) + 2*k2(1) + 2*k3(1) + k4(1));
    y(j+1) = y(j) + (1/6)*h*(k1(2) + 2*k2(2) + 2*k3(2) + k4(2));
    z(j+1) = z(j) + (1/6)*h*(k1(3) + 2*k2(3) + 2*k3(3) + k4(3));
    
    j = j+1;
end

figure(); % Plotting the Lorenz attractor
plot3(x(:,1),y(:,1),z(:,1), 'k',LineWidth=0.6);
xlabel('x'); ylabel('y'); zlabel('z'); view([0 0]);

%% Tangent space
j = 1;
while j<nmax+2
    % k values for perturbation vectors d1 d2 d3
    k1 = lorenz_ts(x(j), y(j), z(j), d1(:,j), d2(:,j), d3(:,j));
    k2 = lorenz_ts(x(j), y(j), z(j), d1(:,j)+h*k1(:,1)/2, d2(:,j)+h*k1(:,2)/2, d3(:,j)+h*k1(:,3)/2);
    k3 = lorenz_ts(x(j), y(j), z(j), d1(:,j)+h*k2(:,1)/2, d2(:,j)+h*k2(:,2)/2, d3(:,j)+h*k2(:,3)/2);
    k4 = lorenz_ts(x(j), y(j), z(j), d1(:,j)+h*k3(:,1), d2(:,j)+h*k3(:,2), d3(:,j)+h*k3(:,3));
    % rk4 applied to perturbation vectors
    A = [d1(:,j) d2(:,j) d3(:,j)] + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);

    % re-orthonormalization
    if rem(j,nnorm) == 0
        [Q, R] = qr(A);
        lam1inst = [lam1inst log(abs(R(1,1)))];
        lam2inst = [lam2inst log(abs(R(2,2)))];
        lam3inst = [lam3inst log(abs(R(3,3)))];
        A = Q;
    end
    d1(:,j+1) = A(:,1); d2(:,j+1) = A(:,2); d3(:, j+1) = A(:,3);
%     mag1(j) = norm(d1(:,j)); mag2(j) = norm(d2(:,j)); mag3(j) = norm(d3(:,j));

    j = j+1;
end

%% Plots

% figure(); plot(mag1); figure(); plot(log(mag1));
% figure(); plot(mag2); figure(); plot(log(mag2));
% figure(); plot(mag3); figure(); plot(log(mag3));

figure(2); %plotting lambda_1 inst wrt time
plot(lam1inst);
ylabel('\lambda_{1i} '); xlabel("time")
figure(3);
plot(cumsum(lam1inst)./(1:length(lam1inst)));
ylabel('Finite-time \lambda_1, \lambda_2, \lambda_3'); xlabel('time');
hold on;
fprintf('Finite time 1st Lyapunov exponent is =')
disp(sum(lam1inst)/(length(lam1inst)));

figure(4); %plotting lambda_2 inst wrt time
plot(lam2inst);
ylabel('\lambda_{2i}'); xlabel("time")
figure(3);
plot(cumsum(lam2inst)./(1:length(lam2inst)));
% ylabel('Finite-time \lambda_2'); xlabel('time');
hold on;
fprintf('Finite time 2nd Lyapunov exponent is')
disp(sum(lam2inst)/(length(lam2inst)));

figure(5); %plotting lambda_3 inst wrt time
plot(lam3inst);
ylabel('\lambda_{3i}'); xlabel("time")
figure(3);
plot(cumsum(lam3inst)./(1:length(lam3inst)));
legend('\lambda_1', '\lambda_2', '\lambda_3');
hold off;
% ylabel('Finite-time \lambda_3'); xlabel('time');
fprintf('Finite time 3rd Lyapunov exponent is')
disp(sum(lam3inst)/length(lam3inst));

%% Function to simulate lorenz
function vec = lorenz_rhs(X, Y, Z)
    sig = 10; r = 28; b = 8/3;
    ex = sig*(Y-X);
    yi = r*X-Y-X*Z;
    zi = X*Y-b*Z;
    vec = [ex yi zi];
end

%% Function for tangent space equations
function vect = lorenz_ts(xlor, ylor, zlor, dX,dY,dZ)
    sig = 10; r = 28; b = 8/3;
    vect = [-sig sig 0; (r-zlor) -1 -xlor; ylor xlor -b]*[dX dY dZ];
end
