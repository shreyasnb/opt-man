function e = checkgradient(problem)
    n = problem.dim1;
    x = rand(n,1);
    t = logspace(-8,0,51);
    v = problem.grad(x);
    f_tv = zeros(length(n),1);
    f_x = zeros(length(n),1);
    e = zeros(length(n),1);
    for i=1:length(t)
        f_tv(i) = problem.cost(x-t(i)*v);
        f_x(i) = problem.cost(x);
        e(i) = abs(f_tv(i) - f_x(i) - t(i).*(v'*v));
    end
    figure
    loglog(t, e, 'Color',[0.5 0.5 0.5]);
    hold on
    grid on
    loglog(t, 2*t, 'LineStyle','--', 'Color','k')
    xlabel('t')
    ylabel('E(t)')
    title('Checkgradient: The grey line should match the slope of the dashed line for a few intervals')
end