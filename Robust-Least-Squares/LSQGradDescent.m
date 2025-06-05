oldpath = path;
path('../utils',oldpath)

setup();
% rng(42, 'twister');

gradf_store = [];
x_store = [];
proj_x_store = [];
f_store = [];
alpha_store = [];

%% Problem initialization
problem.dim = [2,2];
problem.data.theta_0 = 0*pi/180;
problem.data.theta_b = 30*pi/180;
problem.params.alpha = 1e-2;
problem.params.tolx = 1e-12;
problem.params.rho = sin(45*pi/180);
problem.params.lambda = 0.2; 

problem = init(problem);

problem.cost = @(x,y) f(x,problem.Y,problem.data.b);
problem.grad = @(x,y) gradf(x,problem.Y, problem.data.b);
problem.Bx = @(x) B(x, problem);

% checkgradient(problem);

%% Gradient Descent
results = GDA(problem);

%% Plot results
plot_gradf(results);
plot_cost(results, problem);
plot_optimal_sol(results, problem);


%% Functions
function fun = f(x, Y,b)
    fun = norm(Y*(Y'*x) - b)^2;
end

function g = gradf(x, Y, b)
    g = 2*Y*(Y'*(x-b));
end

function Bx = B(x, problem)
    b = problem.data.b;
    Y_hat = problem.data.Y_hat;
    lambda = problem.params.lambda;
    Bx= x*x' - b*x' - x*b' + lambda*(Y_hat*Y_hat');
end

function out = GDA(problem)
    N = 1;
    alpha = problem.params.alpha;
    tolx = problem.params.tolx;
    k = problem.dim(2);
    Y_hat = problem.data.Y_hat;
    b = problem.data.b;
    lambda = problem.params.lambda;
    x_old = problem.x0;

    out.gradf_store = [];
    out.f_store = [];
    out.x_store = [];
    out.proj_x_store = [];

    while (N<=1e5)
        B = problem.Bx(x_old);
        [Y,~] = eigs(B, k);
        v = gradf(x_old, Y,b);
    
        out.gradf_store = [out.gradf_store, v];
        out.f_store = [out.f_store, f(x_old, Y,b)];
        out.x_store = [out.x_store, x_old];
        out.proj_x_store = [out.proj_x_store, Y*(Y'*x_old)];
        
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

    if(N<=1e5)
        disp('Gradient descent converged successfully')
        disp(strcat('Number of iterations:  ',int2str(N)));
        disp(strcat('Norm of gradient:  ', num2str(norm(v))));
    else
        if(norm(v)>tolx)
            disp('Number of iterations exceeded! \n')
            disp('Gradient descent failed within 1e5 iterations');
        end
    end
end