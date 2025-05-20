function fig_attractor(dynamics)
%FIG_ATTRACTOR Summary of this function goes here
%   Detailed explanation goes here
figure(); % Plotting the Lorenz attractor
plot3(dynamics(1,:),dynamics(2,:),dynamics(3,:), 'g',LineWidth=0.3);
xlabel('$x$', 'Interpreter', 'latex'); ylabel('$y$', 'Interpreter', 'latex');
zlabel('$z$', 'Interpreter', 'latex'), %view([0 0]);
set(gca,'TickLabelInterpreter','latex');
end

