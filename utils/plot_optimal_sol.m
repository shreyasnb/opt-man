function p = plot_optimal_sol(x_new, Y, problem)
    theta = problem.params.theta;
    b = problem.data.b;
    theta_0 = problem.data.theta_0;
    Y_hat = problem.data.Y_hat;

    xv = -3:3;
    xs = -1:1;
    figure
    plot(Y_hat(1)*xv, Y_hat(2)*xv,'Color','blue')
    hold on
    grid on
    plot(Y(1)*xv, Y(2)*xv, 'Color','k')
    X_fill = [xs, fliplr(xs)];
    Y_fill = [xs*tan(theta), fliplr(xs*tan(-theta))];
    [X_fill_r, Y_fill_r] = rotate(X_fill, Y_fill, theta_0);
    fill(X_fill_r, Y_fill_r, 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    scatter(b(1), b(2), 80,'blue', 'filled')
    scatter(x_new(1), x_new(2), 60,'k', 'filled');
    [x1,y1] = rotate(Y_hat(1)*xv, Y_hat(2)*xv, theta);
    [x2, y2] = rotate(Y_hat(1)*xv, Y_hat(2)*xv, -theta);
    plot(x1, y1, 'LineStyle','--', 'Color','red')
    plot(x2, y2, 'LineStyle','--', 'Color','red')
    legend({'$\hat{y}$', '$y^*$', '$d_2(y,\hat{y}) = \rho$', '$b$', '$x^*$'},'Interpreter','latex');
    xlim([-1.5 1.5])
    % title('Optimal solution')
end