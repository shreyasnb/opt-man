oldpath = path;
path('../utils',oldpath)

setup();

n = 3;
N = 1;
p=50;

A = randn(n);
A = (A + A')/2;

alpha_0 = 1e-2;
tau = 0.5;
r = 1e-4;
norm_tol = 1e-10;
x_old = [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];

v = f_x(x_old, A);
x_new = retr(x_old, alpha_0*proj_x(x_old,v));

f_old = f(x_old, A);
f_store = [];
gradf_store = [];
alpha_store = [];
alpha = alpha_0;

generate_sphere(p);
hold on

while (N<1e4)
    f_store = [f_store, f(x_old, A)];
    gradf_store = [gradf_store, gradf(x_old, A)];
    alpha_store = [alpha_store, alpha];
    v = gradf(x_old, A);
    alpha = alpha_0;
    while true
        if((f(x_old,A)-f(retr(x_old, alpha*v), A))< r*alpha*norm(v)^2)
            break
        else
            alpha = tau*alpha;
        end
    end
    if ((norm(gradf(x_old, A)) < norm_tol))
        break
    end
    x_old = x_new;
    x_new = retr(x_old, alpha*v); 
    [az,el,r] = cart2sph(x_new(1), x_new(2), x_new(3));
    [xs,ys,zs] = sph2cart(az, el, linspace(0,1,100));
    s = scatter3(xs(end),ys(end),zs(end), 7, "k", "filled", "o");
    set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    N = N + 1;
end

plot_eigvectors(A)

plot_convergence(A, f_store);
plot_gradf(gradf_store);
plot(1:N, alpha_store);

function cost = f(x,A)
    cost = x'*A*x;
end

function egradf = f_x(x,A)
    egradf = 2*A*x;
end

function g = gradf(x,A)
    g = proj_x(x,f_x(x,A));
end