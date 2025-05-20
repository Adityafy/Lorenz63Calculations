%% CLEs version 4

%{
The second approach:
Construct CLVCs at all of the time steps using the stored values of the c matrix at those time steps.
Normalize these CLV’s and you can now make a movie of the CLV’s that is accurate for long times.
Also you can take each one of these normalized CLV's and integrate it forward in time using the Jacobian 
for just one time step and then compute the instantaneous CLE from the magnitude change after that one time step. 
You can then average all of the those instantaneous values to get the finite time CLE.
%}

close all;
clear all;

%% Lorenz solver with RK4

% parameters must be changed in functions if needed
sigma = 10; r = 28; b = 8/3; %parameters

nmax = 200000; %number of time-steps
deltat = 0.001; %time-step length
h = deltat;

ic = [-7.1778; -12.9972; 12.5330]; % initial conditions
X = [ic(1); ic(2); ic(3)];
for j = 1:nmax
    % k values, rk4 for lorenz
    k1 = lorenzDynamics(X(:,j));
    k2 = lorenzDynamics(X(:,j) + (0.5*h)*k1);
    k3 = lorenzDynamics(X(:,j) + (0.5*h)*k2);
    k4 = lorenzDynamics(X(:,j) + h*k3);
    % rk4 applied to lorenz system
    X(:,j+1) = X(:,j) + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
end

% splitting x into x,y,z and making them have exactly nmax values
x = X(1,1:nmax); y = X(2,1:nmax); z = X(3,1:nmax);

figure(); % Plotting the Lorenz attractor
plot3(x(1,:),y(1,:),z(1,:), 'k',LineWidth=0.3);
xlabel('$x$', 'Interpreter', 'latex'); ylabel('$y$', 'Interpreter', 'latex');
zlabel('$z$', 'Interpreter', 'latex'), %view([0 0]);
set(gca,'TickLabelInterpreter','latex');


%% Tangent space, QR decomposition

% d1 = [deltax1; deltay1; deltaz1] (similar for d2, d3)
ic1 = [0; 0; 1]; ic2 = [0; 1; 0]; ic3 = [1; 0; 0]; % initial conditions
d1 = [ic1 zeros(3, nmax)];
d2 = [ic2 zeros(3, nmax)];
d3 = [ic3 zeros(3, nmax)];

nnorm = 1000; % iteration gap for normalization for GS vectors
tnorm = nnorm*deltat;

% initialization of Q and R multidimenstional matrices
Q(:,:,1) = zeros(3,3);
R(:,:,1) = zeros(3,3);

laminst = [];

j = 1;
while j<nmax+1
    
    % k values for perturbation vectors d1 d2 d3
    k1 = lorenz_jacob(x(j),y(j),z(j), d1(:,j), d2(:,j), d3(:,j));
    k2 = lorenz_jacob(x(j),y(j),z(j), d1(:,j)+h*k1(:,1)/2, d2(:,j)+h*k1(:,2)/2, d3(:,j)+h*k1(:,3)/2);
    k3 = lorenz_jacob(x(j),y(j),z(j), d1(:,j)+h*k2(:,1)/2, d2(:,j)+h*k2(:,2)/2, d3(:,j)+h*k2(:,3)/2);
    k4 = lorenz_jacob(x(j),y(j),z(j), d1(:,j)+h*k3(:,1), d2(:,j)+h*k3(:,2), d3(:,j)+h*k3(:,3));
    
    % rk4 applied to perturbation vectors
    A = [d1(:,j) d2(:,j) d3(:,j)] + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);

    % qr decomposition, every time step
    [Qj, Rj] = qr(A);
    Q(:,:,j+1) = Qj;
    R(:,:,j+1) = Rj;

    A = Qj;
    d1(:,j+1) = A(:,1);
    d2(:,j+1) = A(:,2);
    d3(:,j+1) = A(:,3);

    j = j+1;
end

%% checking the rank of R
for j = 1:length(R)
    if rank(R(:,:,j))<3
        disp(j);
    end
end

%% C matrix
nnorm2 = 10;    % time-steps for normalization for clv calculations of c matrix
tnorm2 = nnorm2*deltat;

rng(1);
c1(:,nmax+1) = rand(1,1);
c2(:,nmax+1) = rand(2,1);
c3(:,nmax+1) = rand(3,1);

cmag = zeros(3,1);

for j = nmax:-1:1
    c1(:,j) = inv(R(1,1,j+1))*c1(1,j+1);
    c2(:,j) = inv(R(1:2,1:2,j+1))*c2(1:2,j+1);
    c3(:,j) = inv(R(1:3,1:3,j+1))*c3(1:3,j+1);
    
    cmag(:,j) = [norm(c1(1,j)); norm(c2(:,j)); norm(c3(:,j))];

    if j ~= nmax && rem(j,nnorm2) == 0
        c1(j) = c1(j)/norm(c1(j));
        c2(:,j) = c2(:,j)/norm(c2(:,j));
        c3(:,j) = c3(:,j)/norm(c3(:,j));
    end
end


%% CLVs

nnorm3 = 1; % for clv exponents calculations
tnorm3 = nnorm3*deltat;

w1(:,1) = zeros(3,1);
w2(:,1) = zeros(3,1);
w3(:,1) = zeros(3,1);

% wmag(:,1) = zeros(3,1);
%createw = nmax; % round(nmax/3); % time where we start growing the vectors with Jacobian

for j = 1:nmax
    w1(:,j) = c1(1,j).*d1(:,j);
    w2(:,j) = c2(1,j).*d1(:,j) + c2(2,j).*d2(:,j);
    w3(:,j) = c3(1,j).*d1(:,j) + c3(2,j).*d2(:,j) + c3(3,j).*d3(:,j);
%     wmag(:,j) = [norm(w1(:,j)); norm(w2(:,j)); norm(w3(:,j))];
    if rem(j,nnorm3) == 0
        w1(:,j) = w1(:,j)/norm(w1(:,j));
        w2(:,j) = w2(:,j)/norm(w2(:,j));
        w3(:,j) = w3(:,j)/norm(w3(:,j));
    end
end

% figure(); hold on;
% plot(wmag(1,:)); plot(wmag(2,:)); plot(wmag(3,:));

Wmag(:,1) = zeros(3,1);
cle(:,1) = zeros(3,1);
W(:,:,1) = [w1(:,j) w2(:,j) w3(:,j)];

for j = 1:nmax %createw:nmax

    % k values for perturbation vectors d1 d2 d3
    k1 = lorenz_jacob(x(j),y(j),z(j), w1(:,j), w2(:,j), w3(:,j));
    k2 = lorenz_jacob(x(j),y(j),z(j), w1(:,j)+h*k1(:,1)/2, w2(:,j)+h*k1(:,2)/2, w3(:,j)+h*k1(:,3)/2);
    k3 = lorenz_jacob(x(j),y(j),z(j), w1(:,j)+h*k2(:,1)/2, w2(:,j)+h*k2(:,2)/2, w3(:,j)+h*k2(:,3)/2);
    k4 = lorenz_jacob(x(j),y(j),z(j), w1(:,j)+h*k3(:,1), w2(:,j)+h*k3(:,2), w3(:,j)+h*k3(:,3));
    
    W(:,:,j) = [w1(:,j) w2(:,j) w3(:,j)] + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4);
    
%     Wmag(:,j) = [norm(W(:,1,j))-norm(w1(:,j)); norm(W(:,2,j))-norm(w2(:,j)); norm(W(:,3,j))-norm(w3(:,j))];
    Wmag(:,j) = [norm(W(:,1,j)); norm(W(:,2,j)); norm(W(:,3,j))];
 
%     if rem(j-createw,nnorm3) == 0
%     cle(:,j) = [log(abs(Wmag(1,j))); log(abs(Wmag(2,j))); log(abs(Wmag(3,j)))];
    cle(:,j) = (1/tnorm3)*[log(abs(Wmag(1,j))); log(abs(Wmag(2,j))); log(abs(Wmag(3,j)))];
%         for m = 1:3
%             W(:,m) = W(:,m)/norm(W(:,m));
%         end
%     end
    
%     w1(:,j+1) = W(:,1); w2(:,j+1) = W(:,2); w3(:,j+1) = W(:,3);
end

% cle = cle(:,2:end); % first CLEs are initialized as zero, so ignoring that

%% Plotting finite-time CLEs
figure();
hold on;
%timeaxix = linspace(createw*deltat,nmax*deltat,deltat*(nmax-createw)/nnorm3);
plot(cumsum(cle(1,:))./(1:length(cle(1,:))), 'LineWidth',1);
set(gca,'TickLabelInterpreter','latex');
ylabel('Finite-time CLE $\lambda_1$ , $\lambda_2$ , $\lambda_3$', 'Interpreter','latex', 'FontWeight','bold');
xlabel('Time-steps','Interpreter','latex', 'FontWeight','bold');
plot(cumsum(cle(2,:))./(1:length(cle(2,:))), 'LineWidth',1);
plot(cumsum(cle(3,:))./(1:length(cle(3,:))), 'LineWidth',1);
legend('\lambda_1','\lambda_2','\lambda_3');
% xlim([0 length(cle(1,:))]);
% xticks([createw*deltat deltat*(nmax-createw)/length(cle(1,:)) nmax*deltat]);
hold off;


% %% Plotting the vectors (animation)
% figure();
% set(gca,'TickLabelInterpreter','latex');
% plot3(x(1,:),y(1,:),z(1,:), 'k',LineWidth=0.2);
% xlabel('$x$', 'Interpreter', 'latex');
% ylabel('$y$', 'Interpreter', 'latex');
% zlabel('$z$', 'Interpreter', 'latex'),
% view([0 0]);
% hold on
% scale=10; % scaling for plots
% gsvscale=0.7*scale;
% lw = 4; % linewidth for plots
% headsize = 1.5; %arrowhead size
% 
% for j = 1:length(w1)
% 
%     h1 = plot3(x(j),y(j),z(j),'r.','Markersize', 30);
%     v = lorenzDynamics(X(:,j)); % velocity
%     v = scale*v/norm(v);
%     h2 = quiver3(x(j),y(j),z(j),v(1),v(2),v(3),'k', 'LineWidth',lw,MaxHeadSize=headsize);
% 
%     d1j = gsvscale*(d1(:,j)/norm(d1(:,j)));
%     h3 = quiver3(x(j),y(j),z(j),d1j(1),d1j(2),d1j(3), 'c', 'LineWidth',1.1*lw,MaxHeadSize=headsize);
%     d2j = gsvscale*(d2(:,j)/norm(d2(:,j)));
%     h4 = quiver3(x(j),y(j),z(j),b*d2j(1),b*d2j(2),b*d2j(3), 'm', 'LineWidth',lw,MaxHeadSize=headsize);
%     d3j = gsvscale*(d3(:,j)/norm(d3(:,j)));
%     h5 = quiver3(x(j),y(j),z(j),b*d3j(1),b*d3j(2),b*d3j(3), 'y', 'LineWidth',lw,MaxHeadSize=headsize);
%     
%     w1j = scale*(w1(:,j)/norm(w1(:,j)));
%     w2j = scale*(w2(:,j)/norm(w2(:,j)));
%     w3j = scale*(w3(:,j)/norm(w3(:,j)));
%     h6 = quiver3(x(j),y(j),z(j),w1j(1),w1j(2),w1j(3),'r', 'LineWidth',lw,MaxHeadSize=headsize);
%     h7 = quiver3(x(j),y(j),z(j),w2j(1),w2j(2),w2j(3),'b', 'LineWidth',lw,MaxHeadSize=headsize);
%     h8 = quiver3(x(j),y(j),z(j),w3j(1),w3j(2),w3j(3),'g', 'LineWidth',lw,MaxHeadSize=headsize);
% 
%     drawnow;
% 
%     delete(h1);
%     delete(h2);
%     delete(h3);
%     delete(h4);
%     delete(h5);
%     delete(h6);
%     delete(h7);
%     delete(h8);
% end

%% Plotting the vectors (animation) clv3
% Create figure
vectorPlot = figure;
% Create axes
axes1 = axes('Parent',vectorPlot);
hold(axes1,'on');
startTime = 1200;
endTime = 1600;
plot3(x(1,startTime:endTime+startTime),...
    y(1,startTime:endTime+startTime),...
    z(1,startTime:endTime+startTime),'k',LineWidth=0.2);
% set(gca,'TickLabelInterpreter','latex');
% xlabel('$x$', 'Interpreter', 'latex');
% ylabel('$y$', 'Interpreter', 'latex');
% zlabel('$z$', 'Interpreter', 'latex'),
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
view([30 15]);
xlim([-20 20]);
ylim([-25 25]);
zlim([0 50]);
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',15);
hold on
scale=12; % scaling for plots
gsvscale=0.7*scale;
lw = 4; % linewidth for plots
headsize = 1.5; %arrowhead size


for j = startTime:endTime

    h1 = plot3(x(j),y(j),z(j),'r.','Markersize', 30);
    trajVec = lorenzDynamics(X(:,j)); % velocity
    trajVec = gsvscale*trajVec/norm(trajVec);
    h2 = quiver3(x(j),y(j),z(j),trajVec(1),trajVec(2),trajVec(3), 'LineWidth',lw,MaxHeadSize=headsize,Color=[0.4392 0.7608 0.0235]);

    d1j = gsvscale*(d1(:,j)/norm(d1(:,j)));
    d2j = gsvscale*(d2(:,j)/norm(d2(:,j)));
    d3j = gsvscale*(d3(:,j)/norm(d3(:,j)));

    h3 = quiver3(x(j),y(j),z(j),d1j(1),d1j(2),d1j(3), 'LineWidth',1.4*lw,MaxHeadSize=headsize, Color=[0 0.7882 0.7882]);

    h4 = quiver3(x(j),y(j),z(j),b*d2j(1),b*d2j(2),b*d2j(3), 'LineWidth',lw,MaxHeadSize=headsize,Color=[0.8588 0.4627 0.7020]);

    h5 = quiver3(x(j),y(j),z(j),b*d3j(1),b*d3j(2),b*d3j(3),'LineWidth',lw,MaxHeadSize=headsize,Color=[0.6020 0.6020 0.6020]);
    
    w1j = scale*(w1(:,j)/norm(w1(:,j)));
    w2j = scale*(w2(:,j)/norm(w2(:,j)));
    w3j = scale*(w3(:,j)/norm(w3(:,j)));
    h6 = quiver3(x(j),y(j),z(j),w1j(1),w1j(2),w1j(3), 'LineWidth',lw,MaxHeadSize=headsize,Color=[0.6902 0.0196 0.0196]);
    h7 = quiver3(x(j),y(j),z(j),w2j(1),w2j(2),w2j(3), 'LineWidth',lw,MaxHeadSize=headsize,Color=[0 0.2667 1.0000]);
    h8 = quiver3(x(j),y(j),z(j),w3j(1),w3j(2),w3j(3), 'LineWidth',lw,MaxHeadSize=headsize,Color=[0.8902 0.4980 0.3020]);

    drawnow;
    if j < endTime
        delete(h1);
        delete(h2);
        delete(h3);
        delete(h4);
        delete(h5);
        delete(h6);
        delete(h7);
        delete(h8);
    end
end

%% Function for Lorenz Equations
function vec = lorenzDynamics(X)
    sig = 10; r = 28; b = 8/3;
    x = sig*(X(2,:)-X(1,:));
    y = r*X(1,:)-X(2,:)-X(1,:)*X(3,:);
    z = X(1,:)*X(2,:)-b*X(3,:);
    vec = [x; y; z];
end

%% Jacobian Function for tangent space equations
function mat = lorenz_jacob(xlor, ylor, zlor, dX,dY,dZ)
    sig = 10; r = 28; b = 8/3;
    mat = [-sig sig 0; (r-zlor) -1 -xlor; ylor xlor -b]*[dX dY dZ];
end
