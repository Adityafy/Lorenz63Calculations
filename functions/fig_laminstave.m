function fig_laminstave(laminst)
% fig_laminstave(laminst)
% use any instantaneous \lambda, either from gs or clvs
figure;
hold on;
plot(cumsum(laminst(1,:))./(1:length(laminst(1,:))), ...
    'DisplayName','$\lambda_1$','LineWidth',1.5);
plot(cumsum(laminst(2,:))./(1:length(laminst(2,:))), ...
    'DisplayName','$\lambda_2$','LineWidth',1.5);
plot(cumsum(laminst(3,:))./(1:length(laminst(3,:))), ...
    'DisplayName','$\lambda_3$','LineWidth',1.5);
legend;
xlim([1,length(laminst(1,:))]);
set(legend,'Interpreter','latex');
set(gca,'FontSize',15,'Box','on','XMinorTick','on','YMinorTick','on');
xlabel('$n$','Interpreter','latex','FontSize',30);
ylabel('$\left<\lambda(t)\right>_t$','Interpreter','latex', ...
    'FontSize',30,'Rotation',0);
axis square;
end

