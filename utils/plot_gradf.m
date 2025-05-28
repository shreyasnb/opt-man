function f = plot_gradf(gfs)
    N = length(gfs);
    norm_gfs = zeros(1,N);
    for i=1:N
        norm_gfs(1,i) = norm(gfs(:,i));
    end
    figure
    loglog(1:N, norm_gfs, 'Color', 'blue');
    grid on
    % hold on
    % plot(1:N, e*ones(N), 'Color','red');
    xlabel('Iteration');
    ylabel('grad$(f)$')
    title('Gradient of cost function');
    hold off
    % legend({'$f(x) = x^T A x$', '$\lambda_{max}(A)$'},'Location','best', 'Interpreter','latex')
end