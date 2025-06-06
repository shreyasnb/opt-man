function p = plot_optimal_sol(results, problem)
    theta = asin(problem.params.rho);
    b = problem.data.b;
    theta_0 = problem.data.theta_0;
    Y_hat = problem.data.Y_hat;
    n = problem.dim(1);
    k = problem.dim(2);

    x_store = results.x_store;
    proj_x_store = results.proj_x_store;
    Y_star = results.Y_star;
    lambdas = results.lambda_store;

    x_opt = x_store(:,end);
    proj_x_opt = proj_x_store(:,end);

    figure
    plot(1:length(lambdas), lambdas)
    xlabel('Iteration')
    ylabel('$\lambda_k$')
    grid on
    
    if(n==2 && k==1)
        xv = -2:2;
        xs = -1:1;
        figure
        plot(Y_hat(1)*xv, Y_hat(2)*xv,'Color','blue')
        hold on
        grid on
        plot(Y_star(1)*xv, Y_star(2)*xv, 'Color','k')
        X_fill = [xs, fliplr(xs)];
        Y_fill = [xs*tan(theta), fliplr(xs*tan(-theta))];
        [X_fill_r, Y_fill_r] = rotate(X_fill, Y_fill, theta_0);
        fill(X_fill_r, Y_fill_r, 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        scatter(b(1), b(2), 100,'blue', 'filled')
        N = length(x_store);
        N_inst = 10;
        % for i=0:N_inst
        %     scatter(proj_x_store(1,1+int32(i*(N-1)/N_inst)), proj_x_store(2,1+int32(i*(N-1)/N_inst)), 40, "green", 'filled');
        %     scatter(x_store(1,1+int32(i*(N-1)/N_inst)), x_store(2,1+int32(i*(N-1)/N_inst)), 40, 'yellow','filled','Marker','square');
        %     % alpha(s, tanh((i+1)/(length(x_store)+1)))
        % end
        s=scatter(proj_x_store(1,:), proj_x_store(2,:), 40, "green", 'filled','^');
        % alpha(s, tanh((1:N)/(N+1)));
        scatter(proj_x_opt(1), proj_x_opt(2), 100,'k', 'filled');
        [x1,y1] = rotate(Y_hat(1)*xv, Y_hat(2)*xv, theta);
        [x2, y2] = rotate(Y_hat(1)*xv, Y_hat(2)*xv, -theta);
        plot(x1, y1, 'LineStyle','--', 'Color','red')
        plot(x2, y2, 'LineStyle','--', 'Color','red')
        legend({'$\hat{y}$', '$y^*$','$d_2(y,\hat{y}) = \rho$', '$b$', '$P_{y_k}x_k$', '$P_{y^*}x^*$'},'Interpreter','latex');
        % title('Optimal solution')
        xlim([-1.2 1.2])
        ylim([-1.2 1.2])
    elseif(n==3 && k==1)
        xv = -3:3;
        xs = -1:1;
        figure
        [X,Y,Z] = cylinder([0 tan(theta)], 50);
        axis([-3 3,-3 3,-3 3])
        M=makehgtform('translate',[0,0,0],'xrotate',pi/2,'yrotate',pi/2);
        Mp = makehgtform('translate', [0,0,0], 'xrotate', pi/2, 'yrotate',-pi/2);
        h=surf(X,Y,Z,'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.4);
        hold on
        hp = surf(X,Y,Z, 'Parent',hgtransform('Matrix', Mp), 'LineStyle','none', 'FaceAlpha',0.4);
        view([30,35])
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        grid on
        light
        plot3(Y_hat(1)*xv, Y_hat(2)*xv, Y_hat(3)*xv, 'Color', 'blue');
        hold on
        plot3(Y_star(1)*xv, Y_star(2)*xv, Y_star(3)*xv, 'Color','k');
        scatter3(b(1), b(2),b(3), 80, 'blue', 'filled')
        scatter3(proj_x_opt(1), proj_x_opt(2), proj_x_opt(3), 60, 'k', 'filled');
        legend({'$\hat{y}$','$y^*$','$b$','$x^*$'},'Interpreter','latex');
    end
end