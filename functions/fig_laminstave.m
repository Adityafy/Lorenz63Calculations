function fig_laminstave(laminst)
%FIG_LAMINSTAVE Summary of this function goes here
%   Detailed explanation goes here
figure;
hold on;
plot(cumsum(laminst(1,:))./(1:length(laminst(1,:))), ...
    'DisplayName','$\lambda_1$');
plot(cumsum(laminst(2,:))./(1:length(laminst(2,:))), ...
    'DisplayName','$\lambda_2$');
plot(cumsum(laminst(3,:))./(1:length(laminst(3,:))), ...
    'DisplayName','$\lambda_3$');
legend;
set(legend,'Interpreter','latex');
end

