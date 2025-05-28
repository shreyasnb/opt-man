function s = plot_convergence(A, fs)
    e = max(eigs(A));
    N = length(fs);
    figure
    plot(1:N, fs, 'Color', 'blue');
    grid on
    hold on
    plot(1:N, e*ones(N), 'Color','red');
    xlabel('Iteration');
    ylabel('Maximum Eigenvalue of $A$')
    title('Convergence of optimal cost to maximum eigenvalue');
    legend({'$f(x) = x^T A x$', '$\lambda_{max}(A)$'},'Location','best', 'Interpreter','latex')
end