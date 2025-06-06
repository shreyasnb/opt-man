oldpath = path;
path('../utils',oldpath)

setup();
% rng(42, 'twister');

%% Problem initialization
problem.dim = [2,1];
problem.data.theta_0 = 0*pi/180;
problem.data.theta_b = 25*pi/180;
problem.params.alpha = 1e-2;
problem.params.tolx = 1e-12;
problem.params.rho = sin(20*pi/180);
problem.params.lambda = 0; 

problem = init(problem);

problem.cost = @(x,y) f(x,problem.Y,problem.data.b);
problem.grad = @(x,y) gradf(x,problem.Y, problem.data.b);
problem.Bx = @(x) B(x, problem.params.lambda, problem);

% checkgradient(problem);

%% Gradient Descent
results = GDA(problem);

%% Plot results
plot_gradf(results);
plot_cost(results, problem);
plot_optimal_sol(results, problem);


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

function Bx = B(x, lambda, problem)
    b = problem.data.b;
    Y_hat = problem.data.Y_hat;
    Bx= x*x' - b*x' - x*b' + lambda*(Y_hat*Y_hat');
end

function Ys = top_k(B,k)
    [Ys,~] = eigs(B,k);
end

function ds = subspace(lambda, params)
    rho = params.rho;
    x = params.x;
    Y_hat = params.data.Y_hat;
    k = params.k;

    Bx = B(x,lambda, params);
    ds = d2(top_k(Bx,k),Y_hat) - rho;
end

function out = GDA(problem)
    N = 1;
    alpha = problem.params.alpha;
    tolx = problem.params.tolx;
    rho = problem.params.rho;
    k = problem.dim(2);
    Y_hat = problem.data.Y_hat;
    b = problem.data.b;
    lambda_0 = 0;
    x_old = problem.x0;
    h = 1e-3;

    % params.data.b = b;
    % params.data.Y_hat = Y_hat;
    % params.k = k;
    % params.rho = rho;

    out.gradf_store = [];
    out.f_store = [];
    out.x_store = [];
    out.proj_x_store = [];
    out.lambda_store = [];

    while(N<=1e5)
        lambda_old = lambda_0;
        while true
            Bi = B(x_old, lambda_old, problem);
            [Y,~] = eigs(Bi,k);
            if(d2(Y,Y_hat)<=rho)
                lambda_crit = lambda_old;
                break
            else
                lambda_old = lambda_old + h;
            end
        end
        v = gradf(x_old, Y,b);
    
        out.gradf_store = [out.gradf_store, v];
        out.f_store = [out.f_store, f(x_old, Y, problem)];
        out.x_store = [out.x_store, x_old];
        out.proj_x_store = [out.proj_x_store, Y*(Y'*x_old)];
        out.lambda_store = [out.lambda_store, lambda_crit];
        
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
        disp('Gradient descent converged successfully')
        disp(strcat('Number of iterations:  ',int2str(N)));
        disp(strcat('Norm of gradient:  ', num2str(norm(v))));
    else
        if(norm(v)>tolx)
            disp('Gradient descent failed within 1e5 iterations');
        end
    end
end