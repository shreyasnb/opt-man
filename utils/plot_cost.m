function p = plot_cost(results,problem)
    f_store = results.f_store;
    b = problem.data.b;
    Y = results.Y_star;
    N = length(f_store);
    f_opt = norm(Y*(Y'*b) - b)^2;
    figure
    plot(1:N, f_store, 'DisplayName','$f(x_k,y_k)$');
    hold on
    grid on
    plot(1:N, f_opt*ones(1,N), 'DisplayName','$f^*$');
    title('Cost function')
    xlabel('Iteration')
    ylabel('$f(x,y)$', 'Interpreter','latex');
    legend('Interpreter','latex')
end