oldpath = path;
path('../utils',oldpath)

setup();
% rng(42, 'twister');
n = 2;
k = 1;

tolx = 1e-12;
N = 1;
alpha_0 = 1e-2;
r = 1e-4;
rho = sin(2*pi/8);
theta = asin(rho);
theta_0 = 0;
theta_x = 10*pi/180;
theta_b = 2*pi/16;
lambda = 0.1;
tau = 0.1;

gradf_store = [];
alpha_store = [];
alpha = alpha_0;

problem.dim1 = n;
problem.dim2 = k;
problem.data.theta_0 = theta_0;
problem.data.theta_b = theta_b;
problem.params.rho = rho;
problem.params.theta = theta;

[b, x_old, Y, Y_hat] = init(problem);

problem.cost = @(x) f(x,Y,Y_hat, b, lambda);
problem.grad = @(x) gradf(x,Y,b);
problem.data.b = b;
problem.data.Y_hat = Y_hat;

% checkgradient(problem);

v = gradf(x_old, Y, b);
x_new = x_old - alpha*v;

while (N<=1e5)
    v = gradf(x_old, Y, b);
    gradf_store = [gradf_store, v];
    % alpha_store = [alpha_store, alpha];
    % alpha = alpha_0;
    % while true
    %     if((f(x_old,Y, Y_hat, b, lambda)-f((x_old-alpha*v), Y, Y_hat, b, lambda))< -r*alpha*norm(v)^2)
    %         break
    %     else
    %         alpha = tau*alpha;
    %     end
    % end
    if(norm(v)< tolx)
        break
    else
        x_old = x_new;
        x_new = x_old - alpha*v;
        Bx = B(x_new, b, Y_hat, lambda);
        [Y,~] = eigs(Bx, k);
    end
    N = N+1;
end

plot_gradf(gradf_store);
plot_optimal_sol(x_new, Y, problem);

function fun = f(x, Y, Y_hat, b, lambda)
    fun = trace(Y'*B(x,b,Y_hat, lambda)*Y);
end

function g = gradf(x, Y, b)
    g = 2*(Y*Y')*(x-b);
end

function Bout = B(x, b, Y_hat, lambda)
    Bout = x*x' - b*x' - x*b' + lambda*(Y_hat*Y_hat');
end