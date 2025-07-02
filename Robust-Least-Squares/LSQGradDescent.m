oldpath = path;
path('../utils',oldpath)

setup();
    
%% Problem initialization
problem_description.dim = [2,1];
problem_description.data.theta_0 = 0*pi/180;
problem_description.data.theta_b = 35*pi/180;
problem_description.params.alpha = 1e-2;
problem_description.params.tolx = 1e-12;
problem_description.params.rho = sin(5*pi/180);
problem_description.params.lambda = 0; 

%% Algorithm
% problem_description.algorithm = 'eigm';
% problem_description.algorithm = 'eigf';
problem_description.algorithm = 'svd';

%% Flags
problem_description.plots = false;
problem_description.debug = false;

if (problem_description.debug)

    algos = {'eigm','eigf','svd'};
    
    problem_description.algorithm = algos(1);
    [results,~] = robust_lsq(problem_description);  
    gfs = results.gradf_store;
    N1 = length(gfs);
    norms1 = zeros(1,N1);
    for j=1:N1
        norms1(1,j) = norm(gfs(:,j));
    end
    
    problem_description.algorithm = algos(2);
    [results,~] = robust_lsq(problem_description);  
    gfs = results.gradf_store;
    N2 = length(gfs);
    norms2 = zeros(1,N2);
    for j=1:N2
        norms2(1,j) = norm(gfs(:,j));
    end
    
    problem_description.algorithm = algos(3);
    [results,~] = robust_lsq(problem_description);  
    gfs = results.gradf_store;
    N3 = length(gfs);
    norms3 = zeros(1,N3);
    for j=1:N3
        norms3(1,j) = norm(gfs(:,j));
    end
    
    figure
    grid on
    loglog(1:N1, norms1, 'DisplayName',string(algos(1)))
    hold on
    loglog(1:N2, norms2, 'DisplayName', string(algos(2)))
    loglog(1:N3, norms3, 'DisplayName', string(algos(3)))
    legend('Location','best')
    xlabel('Iteration')
    ylabel('grad$(f)$')
else
    [results,problem] = robust_lsq(problem_description);
end

function [results, problem] = robust_lsq(problem_description)
    
    problem = init(problem_description);
    
    problem.cost = @(x,y) f(x,problem.Y,problem.data.b);
    problem.grad = @(x,y) gradf(x,problem.Y, problem.data.b);
    problem.Bx = @(x) B(x, problem.params.lambda, problem);
    
    % checkgradient(problem);
    
    % Gradient Descent
    % tic
    % results = GDA(problem);
    % toc
    
    %% Plot results
    if(problem.plots)
        plot_gradf(results);
        plot_cost(results, problem);
        plot_optimal_sol(results, problem);
    end

    results = analyze_lambda(problem);
    
    %% Functions
    function fun = f(x, Y, problem)
        b = problem.data.b;
        Y_hat = problem.data.Y_hat;
        lambda = problem.params.lambda;
        fun = norm(Y*(Y'*x) - b)^2 + lambda*trace((Y*Y')*(Y_hat*Y_hat'));
    end
    
    function g = gradf(x, Y, b)
        g = 2*Y*(Y'*(x-b));
    end
    
    function d = d2(Y1,Y2)
        n = size(Y1,1);
        d = sqrt(trace((eye(n) - Y1*Y1')*(Y2*Y2')));
    end
    
    function d = d_infty(Y1, Y2)
        d = norm(Y1*Y1' - Y2*Y2');
    end
    
    function Bx = B(x, lambda, problem)
        b = problem.data.b;
        Y_hat = problem.data.Y_hat;
        Bx= x*x' - b*x' - x*b' + lambda*(Y_hat*Y_hat');
    end
    
    function Ax = A(x,b)
        Ax = x*x' - b*x' - x*b';
    end
    
    function l = lambda(x,Y,b)
        theta_x = atan2(x(2),x(1));
        theta_y = atan2(Y(2),Y(1));
        theta_b = atan2(b(2),b(1));
    
        l = ((cos(theta_y-theta_x)-cos(theta_y-theta_b))*cos(theta_x)*(tan(theta_x)-tan(theta_y)) - cos(theta_y-theta_x)*cos(theta_b)*(tan(theta_b)-tan(theta_y)))/sin(theta_y);
    end
    
    function results = analyze_lambda(problem)
        rho = problem.params.rho;
        n = problem.dim(1);
        k = problem.dim(2);
        Y_hat = problem.data.Y_hat;
        h=1e-3;
        b = problem.data.b;
        lambda_0 = problem.params.lambda;
        
        N=100;
        theta_x = linspace(0,2*pi,N);
        rhos = linspace(5*pi/180,pi/4,5);
        theta_bs = linspace(5*pi/180,pi/4,5);
        lambda_store = [];
        Y_store = [];
        
        for j=1:length(theta_bs)
            b = [cos(theta_bs(j)); sin(theta_bs(j))];
            for i=1:length(theta_x)
                x = [cos(theta_x(i)); sin(theta_x(i))];
                lambda_old = lambda_0;
                while true
                    H = [sqrt(lambda_old)*Y_hat, (x-b), 1i*b];
                    [P,~,~] = svd(H);
                    Y = P(:,1:k);
                    if(d2(Y,Y_hat)<=rho)
                        lambda_crit = lambda_old;
                        Y_crit = Y;
                        break
                    else
                        lambda_old = lambda_old + h;
                    end
                end
                lambda_store = [lambda_store, lambda_crit];
                Y_store = [Y_store, Y_crit];
            end 
        end
        results.theta_bs = theta_bs;
        results.rhos = rhos;
        results.lambda_store = lambda_store;
        results.Y_store = Y_store;
        figure
        for i=1:length(theta_bs)
            Yi = Y_store(:,1+(i-1)*N:N*i);
            theta_y = atan2(Yi(2,:),Yi(1,:));
            hold on
            % plot(theta_x.*180/pi, lambda_store(1+(i-1)*N:N*i).*cos(theta_y)./(cos(theta_x-theta_y)), 'DisplayName',strcat('$\theta_b = $', num2str(theta_bs(i)*180/pi)))
            plot(theta_x.*180/pi, theta_y, 'DisplayName',strcat('$\theta_b$ = ', num2str(theta_bs(i)*180/pi)));
        end
        xlabel('$\theta_x$')
        ylabel('$\theta_y$')
        legend('Interpreter', 'latex', 'Location','best')
    end
    
    function out = GDA(problem)
        N = 1;
        alpha = problem.params.alpha;
        tolx = problem.params.tolx;
        rho = problem.params.rho;
        n = problem.dim(1);
        k = problem.dim(2);
        Y_hat = problem.data.Y_hat;
        b = problem.data.b;
        lambda_0 = problem.params.lambda;
        lambda_crit = lambda_0;
        x_old = problem.x0;
        [Y,~] = eigs(B(x_old, lambda_0, problem), k);
        h = 10;
    
        % params.data.b = b;
        % params.data.Y_hat = Y_hat;
        % params.k = k;
        % params.rho = rho;
    
        out.gradf_store = [];
        out.f_store = [];
        out.x_store = [];
        out.proj_x_store = [];
        out.lambda_store = [];
        out.y_store = [];
    
        while(N<=1e5)
            lambda_old = lambda_0;
            while true
            % lambda_old = lambda(x_old,Y,b);
                if(strcmp(problem.algorithm,'eigf'))
                    Bfun = @(u) (x_old'*u)*x_old - (x_old'*u)*b - (b'*u)*x_old + lambda_old*(Y_hat*(Y_hat'*u));
                    [Y,~] = eigs(Bfun,n,k);
                elseif(strcmp(problem.algorithm, 'eigm'))
                    [Y,~] = eigs(B(x_old, lambda_old, problem), k);
                else
                    H = [sqrt(lambda_old)*Y_hat, (x_old-b), 1i*b];
                    [P,~,~] = svd(H);
                    Y = P(:,1:k);
                end
                if(d2(Y,Y_hat)<=rho)
                    lambda_crit = lambda_old;
                    break
                else
                    lambda_old = lambda_old + h;
                    % break
                end
            end
            v = gradf(x_old, Y,b);
        
            out.gradf_store = [out.gradf_store, v];
            out.f_store = [out.f_store, f(x_old, Y, problem)];
            out.x_store = [out.x_store, x_old];
            out.proj_x_store = [out.proj_x_store, Y*(Y'*x_old)];
            out.lambda_store = [out.lambda_store, lambda_crit];
            out.y_store = [out.y_store, Y];
            if(norm(v)< tolx)
                break
            else
                x_new = x_old - alpha*v;
                x_old = x_new;
            end
            N = N+1;
        end
        out.Y_star = Y;
        out.x_star = x_old;
    
        if(d2(Y,Y_hat) > rho)
                disp('Error: Optimal subspace outside ball')
        end
        if(N<=1e5 && norm(gradf(x_old, Y,b))<tolx)
            disp('Gradient descent converged successfully');
            disp(strcat('Number of iterations:  ',int2str(N)));
            disp(strcat('Norm of gradient:  ', num2str(norm(v))));
        else
            if(norm(v)>tolx)
                disp('Gradient descent failed within 1e5 iterations');
            end
        end
    end
end