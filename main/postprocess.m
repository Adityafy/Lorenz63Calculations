

addpath('../functions/');

%% CLE convergence figure
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

%% Lorenz CLV video

% Create figure
vectorPlot = figure;
% Create axes
axes1 = axes('Parent',vectorPlot);
hold(axes1,'on');
startTime = 1200;
endTime = 1600;
plot3(dynamics(1,startTime:endTime+startTime),...
    dynamics(2,startTime:endTime+startTime),...
    dynamics(3,startTime:endTime+startTime),'k',LineWidth=0.2);
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


for t = startTime:endTime

    h1 = plot3(dynamics(1,t),dynamics(2,t),dynamics(3,t),'r.','Markersize', 30);
    trajVec = lorenz63rhs(p,dynamics(:,t)); % velocity
    trajVec = gsvscale*trajVec/norm(trajVec);
    h2 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        trajVec(1),trajVec(2),trajVec(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.4392 0.7608 0.0235]);

    d1j = gsvscale*(v(:,1,t)/norm(v(:,1,t)));
    d2j = gsvscale*(v(:,2,t)/norm(v(:,2,t)));
    d3j = gsvscale*(v(:,3,t)/norm(v(:,3,t)));

    h3 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        d1j(1),d1j(2),d1j(3), 'LineWidth',1.4*lw,MaxHeadSize=headsize, ...
        Color=[0 0.7882 0.7882]);

    h4 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        b*d2j(1),b*d2j(2),b*d2j(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.8588 0.4627 0.7020]);

    h5 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        b*d3j(1),b*d3j(2),b*d3j(3),'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.6020 0.6020 0.6020]);
    
    w1j = scale*(wf(:,1,t)/norm(wf(:,1,t)));
    w2j = scale*(wf(:,2,t)/norm(wf(:,2,t)));
    w3j = scale*(wf(:,3,t)/norm(wf(:,3,t)));
    h6 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        w1j(1),w1j(2),w1j(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.6902 0.0196 0.0196]);
    h7 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        w2j(1),w2j(2),w2j(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0 0.2667 1.0000]);
    h8 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        w3j(1),w3j(2),w3j(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.8902 0.4980 0.3020]);

    drawnow;
    if t < endTime
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

%% Rossler CLV video
close all;
% Create figure
vectorPlot = figure;
% Create axes
axes1 = axes('Parent',vectorPlot);
hold(axes1,'on');
startTime = 1200;
endTime = 3000;
% plot3(dynamics(1,:),dynamics(2,:),dynamics(3,:),'k');
plot3(dynamics(1,startTime:endTime+3*startTime),...
    dynamics(2,startTime:endTime+3*startTime),...
    dynamics(3,startTime:endTime+3*startTime),'k',LineWidth=0.2);
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


for t = startTime:endTime

    h1 = plot3(dynamics(1,t),dynamics(2,t),dynamics(3,t),'r.','Markersize', 30);
    trajVec = rosslerrhs(p,dynamics(:,t)); % velocity
    trajVec = gsvscale*trajVec/norm(trajVec);
    h2 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        trajVec(1),trajVec(2),trajVec(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.4392 0.7608 0.0235]);

    d1j = gsvscale*(v(:,1,t)/norm(v(:,1,t)));
    d2j = gsvscale*(v(:,2,t)/norm(v(:,2,t)));
    d3j = gsvscale*(v(:,3,t)/norm(v(:,3,t)));

    h3 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        d1j(1),d1j(2),d1j(3), 'LineWidth',1.4*lw,MaxHeadSize=headsize, ...
        Color=[0 0.7882 0.7882]);

    h4 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        b*d2j(1),b*d2j(2),b*d2j(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.8588 0.4627 0.7020]);

    h5 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        b*d3j(1),b*d3j(2),b*d3j(3),'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.6020 0.6020 0.6020]);
    
    w1j = scale*(wf(:,1,t)/norm(wf(:,1,t)));
    w2j = scale*(wf(:,2,t)/norm(wf(:,2,t)));
    w3j = scale*(wf(:,3,t)/norm(wf(:,3,t)));
    h6 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        w1j(1),w1j(2),w1j(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.6902 0.0196 0.0196]);
    h7 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        w2j(1),w2j(2),w2j(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0 0.2667 1.0000]);
    h8 = quiver3(dynamics(1,t),dynamics(2,t),dynamics(3,t), ...
        w3j(1),w3j(2),w3j(3), 'LineWidth',lw,MaxHeadSize=headsize, ...
        Color=[0.8902 0.4980 0.3020]);

    drawnow;
    if t < endTime
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